"""Reusable subprocess transport with one cancellable request at a time."""
import json
import signal

from PySide6.QtCore import QObject, QProcess, Signal

PREFIX = b"@poreckit-result:"


class JobProcess(QObject):
    readyReadStandardOutput = Signal()
    finished = Signal(int, object)
    errorOccurred = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.child = QProcess(self)
        self.child.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        self.child.started.connect(self.send)
        self.child.readyReadStandardOutput.connect(self.read)
        self.child.finished.connect(self.exited)
        self.child.errorOccurred.connect(self.error)
        self.buffer = bytearray()
        self.output = bytearray()
        self.request = None
        self.environment = None

    def setProcessChannelMode(self, mode):
        self.child.setProcessChannelMode(mode)

    def setProcessEnvironment(self, environment):
        self.environment = environment

    def start(self, program, args):
        if self.request is not None:
            raise RuntimeError("A plotting request is already running")
        self.request = {"job": args[-2], "result": args[-1]}
        if self.child.state() == QProcess.ProcessState.NotRunning:
            self.buffer.clear()
            if self.environment is not None:
                self.child.setProcessEnvironment(self.environment)
            self._prepare_group()
            self.child.start(program, [*args[:-2], "--serve"])
        else:
            self.send()

    def _prepare_group(self):
        """Run the worker as its own process-group leader so cancelling also
        stops the external tools (samtools/bcftools) the worker starts."""
        import os
        if hasattr(self.child, "setChildProcessModifier"):
            try:
                self.child.setChildProcessModifier(os.setsid)
            except (RuntimeError, OSError):
                pass

    def send(self):
        if self.request is not None:
            self.child.write((json.dumps(self.request) + "\n").encode("utf-8"))

    def read(self):
        self.buffer.extend(bytes(self.child.readAllStandardOutput()))
        while b"\n" in self.buffer:
            line, _, rest = self.buffer.partition(b"\n")
            self.buffer = bytearray(rest)
            if line.startswith(PREFIX):
                try:
                    status = json.loads(line[len(PREFIX):])
                except ValueError:
                    self.output.extend(line + b"\n")
                    self.readyReadStandardOutput.emit()
                    continue
                if self.request is not None:
                    self.request = None
                    self.finished.emit(status["code"], QProcess.ExitStatus.NormalExit)
            else:
                self.output.extend(line + b"\n")
                self.readyReadStandardOutput.emit()

    def readAllStandardOutput(self):
        output = bytes(self.output)
        self.output.clear()
        return output

    def exited(self, code, status):
        self.read()
        if self.buffer:
            self.output.extend(self.buffer)
            self.buffer.clear()
            self.readyReadStandardOutput.emit()
        if self.request is not None:
            self.request = None
            self.finished.emit(code if code else -1, status)

    def errorString(self):
        return self.child.errorString()

    def error(self, error):
        if error == QProcess.ProcessError.FailedToStart:
            self.request = None
        self.errorOccurred.emit(error)

    def kill(self):
        self._terminate(graceful=True)

    def _signal_group(self, sig):
        """Signal the worker's whole process group; False when unavailable.
        Refuses to signal our own group in case the worker did not get its
        own (child modifier unavailable)."""
        import os
        pid = self.child.processId()
        if pid <= 0 or os.name != "posix" or not hasattr(os, "killpg"):
            return False
        try:
            group = os.getpgid(pid)
            if group == os.getpgrp():
                return False
            os.killpg(group, sig)
            return True
        except (ProcessLookupError, PermissionError, OSError):
            return False

    def _terminate(self, graceful):
        from PySide6.QtCore import QTimer
        if self._signal_group(signal.SIGTERM if graceful else signal.SIGKILL):
            if graceful:
                # External tools get the worker's finally blocks via SIGTERM;
                # escalate to SIGKILL if it does not finish in time.
                QTimer.singleShot(2000, lambda: self._terminate(graceful=False))
            return
        import os
        if os.name == "nt" and graceful and self.child.processId() > 0:
            import subprocess
            subprocess.run(["taskkill", "/PID", str(self.child.processId()), "/T", "/F"],
                           capture_output=True, creationflags=0x08000000)
            return
        self.child.kill()

    def waitForFinished(self, timeout):
        return self.child.waitForFinished(timeout)

    def processId(self):
        return self.child.processId()

    def shutdown(self):
        self.blockSignals(True)
        self.request = None
        if self.child.state() != QProcess.ProcessState.NotRunning:
            self._terminate(graceful=True)
            self.child.waitForFinished(3000)
            if self.child.state() != QProcess.ProcessState.NotRunning:
                self._terminate(graceful=False)
                self.child.waitForFinished(1000)
