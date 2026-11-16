#ifndef AppVersion
  #error AppVersion must be supplied by the build script
#endif
#ifndef SourceDir
  #error SourceDir must be supplied by the build script
#endif
#ifndef OutputDir
  #error OutputDir must be supplied by the build script
#endif
[Setup]
AppId={{79314D14-CFD0-4EAF-9972-3741F4FE8795}
AppName=poreckit GUI
AppVersion={#AppVersion}
AppPublisher=poreckit
DefaultDirName={localappdata}\Programs\poreckit-gui
DefaultGroupName=poreckit GUI
PrivilegesRequired=lowest
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
OutputDir={#OutputDir}
OutputBaseFilename=poreckit-gui-{#AppVersion}-windows-x86_64-setup
SetupIconFile=poreckit.ico
UninstallDisplayIcon={app}\poreckit-gui.exe
Compression=lzma2
SolidCompression=yes
WizardStyle=modern
DisableProgramGroupPage=yes
CloseApplications=yes
[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; Flags: unchecked
[Files]
Source: "{#SourceDir}\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs
[Icons]
Name: "{group}\poreckit GUI"; Filename: "{app}\poreckit-gui.exe"
Name: "{autodesktop}\poreckit GUI"; Filename: "{app}\poreckit-gui.exe"; Tasks: desktopicon
[Run]
Filename: "{app}\poreckit-gui.exe"; Description: "Launch poreckit GUI"; Flags: nowait postinstall skipifsilent
