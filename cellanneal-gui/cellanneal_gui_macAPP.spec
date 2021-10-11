# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(['cellanneal_gui.py'],
             pathex=['/Users/lbuchauer/owncube/cellanneal/cellanneal-gui'],
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

a.datas += [('logo_orange.png', '/Users/lbuchauer/owncube/cellanneal/cellanneal-gui/logo_orange.png', 'img'),
           ('cellanneal_button.png', '/Users/lbuchauer/owncube/cellanneal/cellanneal-gui/cellanneal_button.png', 'img'),
           ('logo.icns', '/Users/lbuchauer/owncube/cellanneal/cellanneal-gui/logo.icns', 'img')]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='cellanneal_gui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False,
          disable_windowed_traceback=False,
          target_arch=None,
          codesign_identity=None,
          entitlements_file=None )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='cellanneal_gui')
app = BUNDLE(coll,
             name='cellanneal_gui.app',
             icon='logo.icns',
             bundle_identifier=None)
