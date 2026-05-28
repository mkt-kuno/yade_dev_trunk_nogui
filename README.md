# YADE nogui ブランチ

本ブランチは、X11/ディスプレイ環境を持たない HPC クラスタやサーバ向けに YADE をヘッドレス（GUI なし）でビルド・実行できるよう最適化したものです。

---

## 概要

通常の YADE ビルドでは OpenGL・Qt5・GLUT・QGLViewer・GL2PS などの GUI 依存ライブラリが必要ですが、本ブランチで追加した `YADE_HEADLESS` CMake オプションを使うと、これらをすべて無効化した状態でビルドできます。

主な変更点：

| 変更箇所 | 内容 |
|---|---|
| `CMakeLists.txt` | `YADE_HEADLESS` オプションを追加。ON にすると `ENABLE_GUI`、`ENABLE_GL2PS`、`USE_QT5` を自動的に OFF に設定する。 |
| `cMake/YadePythonHelpers.cmake` | `YADE_HEADLESS=ON` 時、GUI 系 Python モジュール（`pygraphviz`、`Xlib`、`tkinter`）の検索をスキップする。 |
| `py/CMakeLists.txt` | `ENABLE_GUI` が OFF のとき Qt5 インクルードディレクトリの設定をスキップする。 |
| `py/plot.py` | ディスプレイが利用できない場合のみ `matplotlib` バックエンドを `Agg` に切り替え、`mtTkinter` の読み込みを回避する。 |

---

## ビルド手順

### 前提

- CMake 3.14 以上
- C++17 対応コンパイラ（GCC 8+ / Clang 8+ など）
- Python 3.x
- 必要な依存ライブラリ（`libeigen3-dev`、`libboost-all-dev`、`python3-numpy`、`python3-matplotlib`、`python3-sphinx` など）

### ヘッドレスビルド（本ブランチ推奨）

```bash
mkdir build && cd build
cmake .. \
    -DYADE_HEADLESS=ON \
    -DCMAKE_INSTALL_PREFIX=/usr/local
make -j$(nproc)
sudo make install
```

`YADE_HEADLESS=ON` を指定することで、以下が自動的に設定されます：

- `ENABLE_GUI=OFF`
- `ENABLE_GL2PS=OFF`
- `USE_QT5=OFF`

X11 や OpenGL、Qt5 などのディスプレイ関連ライブラリは不要です。

### 通常ビルド（GUI あり）

GUI を含む標準ビルドを行う場合は `YADE_HEADLESS` を指定しないでください：

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make -j$(nproc)
sudo make install
```

---

## 実行方法

インストール後、通常と同じコマンドでシミュレーションスクリプトを実行できます：

```bash
yade myscript.py
```

ヘッドレス環境ではディスプレイへのアクセスが行われないため、`matplotlib` は自動的に非インタラクティブバックエンド（`Agg`）を使用します。画像ファイルへの保存は通常通り可能です。

---

## 注意事項

- `YADE_HEADLESS=ON` でビルドした場合、GUI（3D ビュー、インスペクタ等）は使用できません。
- `plot.py` の Tkinter 依存部分はヘッドレスモードでは自動的にスキップされます。
- GUI が必要な場合は本ブランチではなく上流の [yade-dev/trunk](https://gitlab.com/yade-dev/trunk) を使用してください。

---

## 上流プロジェクト

本ブランチは [yade-dev/trunk](https://gitlab.com/yade-dev/trunk) をベースにしています。YADE 全般のドキュメントや情報については上流を参照してください：

- [オンラインドキュメント](https://www.yade-dem.org/doc/)
- [インストール手順](https://yade-dem.org/doc/installation.html)
