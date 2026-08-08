[Setup]
AppName=Noise Repellent
AppVersion=0.3.1
AppPublisher=Luciano Dato
DefaultDirName={commoncf}\VST3
DisableDirPage=yes
OutputBaseFilename=NoiseRepellent-Win64-Installer
Compression=lzma2/ultra64
SolidCompression=yes
ArchitecturesInstallIn64BitMode=x64

[Files]
; VST3 File
Source: "..\..\staging\Noise Repellent.vst3\*"; DestDir: "{commoncf}\VST3\Noise Repellent.vst3"; Flags: ignoreversion recursesubdirs createallsubdirs

; LV2 File (Installed to Common Files\LV2)
Source: "..\..\staging\Noise Repellent.lv2\*"; DestDir: "{commoncf}\LV2\Noise Repellent.lv2"; Flags: ignoreversion recursesubdirs createallsubdirs
