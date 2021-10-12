# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(['cellgui.py'],
             pathex=['C:\\Users\\yotam\\Desktop\\cellanneal-master\\cellanneal-gui'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             hooksconfig={},
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

a.datas += [('logo_orange.png', 'C:\\Users\\yotam\\Desktop\\cellanneal-master\\cellanneal-gui\\logo_orange.png', 'img'),
   ('cellanneal_button.png', 'C:\\Users\\yotam\\Desktop\\cellanneal-master\\cellanneal-gui\\cellanneal_button.png', 'img'),
   ('logo.ico', 'C:\\Users\\yotam\\Desktop\\cellanneal-master\\cellanneal-gui\\logo.ico', 'img')]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='cellgui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None,
          icon="logo.ico")
