# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(['cellgui.py'],
             pathex=['/Users/lisa/X/lisabu/cellanneal_gui_dev/cellanneal/cellanneal-gui'],
             binaries=[],
             datas=[],
             runtime_hooks=['/Users/lisa/X/lisabu/cellanneal_gui_dev/cellanneal/cellanneal-gui/use_lib.py'],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)

a.datas += [('logo_orange.png', '/Users/lisa/X/lisabu/cellanneal_gui_dev/cellanneal/cellanneal-gui/logo_orange.png', 'img'),
           ('cellanneal_button.png', '/Users/lisa/X/lisabu/cellanneal_gui_dev/cellanneal/cellanneal-gui/cellanneal_button.png', 'img'),
           ('logo.icns', '/Users/lisa/X/lisabu/cellanneal_gui_dev/cellanneal/cellanneal-gui/logo.icns', 'img')]

pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)

exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='cellgui',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True,
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
               name='cellgui')
