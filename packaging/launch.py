if __name__ == "__main__":
    import multiprocessing
    import os
    import sys
    multiprocessing.freeze_support()
    # Windowed Windows executables have no Python stdio by default. QProcess
    # provides inherited pipe handles for our persistent rendering worker.
    if sys.platform == "win32" and sys.argv[1:2] == ["--worker"]:
        import ctypes
        import msvcrt
        kernel = ctypes.WinDLL("kernel32", use_last_error=True)
        kernel.GetStdHandle.argtypes = [ctypes.c_uint32]
        kernel.GetStdHandle.restype = ctypes.c_void_p
        kernel.GetCurrentProcess.restype = ctypes.c_void_p
        kernel.DuplicateHandle.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p,
                                          ctypes.POINTER(ctypes.c_void_p), ctypes.c_uint32, ctypes.c_int, ctypes.c_uint32]
        process = kernel.GetCurrentProcess()
        for name, number, mode in [("stdin", -10, "r"), ("stdout", -11, "w"), ("stderr", -12, "w")]:
            if getattr(sys, name) is None:
                handle = kernel.GetStdHandle(number & 0xffffffff)
                if handle and handle != ctypes.c_void_p(-1).value:
                    copied = ctypes.c_void_p()
                    if not kernel.DuplicateHandle(process, handle, process, ctypes.byref(copied), 0, False, 2):
                        raise ctypes.WinError(ctypes.get_last_error())
                    flags = os.O_RDONLY if mode == "r" else os.O_WRONLY
                    fd = msvcrt.open_osfhandle(copied.value, flags | os.O_BINARY)
                    setattr(sys, name, os.fdopen(fd, mode, encoding="utf-8", buffering=1))
                    setattr(sys, "__"+name+"__", getattr(sys, name))
    from poreckit_gui.app import main
    main()
