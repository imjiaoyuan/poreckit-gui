import os
import tempfile

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="poreckit-mpl-test-"))
os.environ.setdefault("XDG_CONFIG_HOME", tempfile.mkdtemp(prefix="poreckit-config-test-"))
