#!/usr/bin/env bash
# 递归解压 生物信息学实验数据.rar 及其中嵌套的压缩包到 dataset/ (macOS / Linux 通用)
#
# 注意: 本目录含 .codebuddy 等隐藏目录时不受影响, 因为脚本 cd 到自身所在目录.
set -euo pipefail

cd "$(dirname "$0")"

ARCHIVE="生物信息学实验数据.rar"
OUTDIR="dataset"

if [ ! -f "$ARCHIVE" ]; then
    echo "错误: 找不到文件 $ARCHIVE" >&2
    exit 1
fi

mkdir -p "$OUTDIR"

# ---------- 工具探测 ----------
have() { command -v "$1" >/dev/null 2>&1; }

RAR_TOOL=""
if   have unrar; then RAR_TOOL="unrar"
elif have 7z;    then RAR_TOOL="7z"
elif have 7zz;   then RAR_TOOL="7zz"
elif have unar;  then RAR_TOOL="unar"
elif have bsdtar; then RAR_TOOL="bsdtar"
fi

BSDTAR_FALLBACK=""
have bsdtar && BSDTAR_FALLBACK="bsdtar"

# 解开归档类压缩包 (rar/zip/7z/tar.*)
extract_archive() {
    local f="$1" dest="$2"
    case "$RAR_TOOL" in
        unrar)  unrar x -o+ "$f" "$dest/" ;;
        7z)     7z x -y -o"$dest" "$f" ;;
        7zz)    7zz x -y -o"$dest" "$f" ;;
        unar)   unar -o "$dest" "$f" ;;
        bsdtar) bsdtar -xf "$f" -C "$dest" ;;
        *)      return 1 ;;
    esac
}

# ---------- 第一层: 解 RAR ----------
echo "==> 解压 $ARCHIVE -> $OUTDIR/"
if [ -z "$RAR_TOOL" ]; then
    echo "错误: 未找到可用解压工具。请安装 unrar / 7z / unar:" >&2
    echo "  macOS: brew install unrar  或  brew install p7zip" >&2
    echo "  Debian/Ubuntu: sudo apt install unrar-free  或  sudo apt install p7zip-full" >&2
    exit 1
else
    echo "使用 $RAR_TOOL"
    extract_archive "$ARCHIVE" "$OUTDIR"
fi

# ---------- 递归解压嵌套压缩包 ----------
# 思路: 反复扫描 OUTDIR, 每次解一层, 直到本轮没有任何可解压文件为止.
#       gz/bz2/xz 为单文件压缩 -> 解到同目录并按去掉扩展名命名, 然后删除原压缩文件.
round=0
while :; do
    round=$((round + 1))
    changed=0

    # 收集待处理文件 (逐层处理, 处理完一层后重新扫描)
    while IFS= read -r -d '' f; do
        dir=$(dirname "$f")
        base=$(basename "$f")

        case "$base" in
            *.gz|*.bz2|*.xz)
                # 跳过 tar.gz 等, 交由归档分支处理
                case "$base" in
                    *.tar.gz|*.tar.bz2|*.tar.xz) continue ;;
                esac
                out="${f%.*}"
                # 已解过则跳过
                if [ -e "$out" ]; then continue; fi
                echo "  [gzip/bzip2/xz] $f"
                case "$base" in
                    *.gz)  gzip -dc  "$f" > "$out" ;;
                    *.bz2) bzip2 -dc "$f" > "$out" ;;
                    *.xz)  xz -dc    "$f" > "$out" ;;
                esac
                rm -f "$f"
                changed=1
                ;;
            *.tar.gz|*.tgz|*.tar.bz2|*.tbz2|*.tar.xz|*.txz|*.tar)
                staged="${f%.tar.gz}"; staged="${staged%.tgz}"
                staged="${staged%.tar.bz2}"; staged="${staged%.tbz2}"
                staged="${staged%.tar.xz}"; staged="${staged%.txz}"
                staged="${staged%.tar}"
                if [ -e "$staged" ]; then continue; fi
                echo "  [tar] $f"
                mkdir -p "$staged"
                bsdtar -xf "$f" -C "$staged"
                rm -f "$f"
                changed=1
                ;;
            *.zip)
                stage="${f%.zip}"
                if [ -e "$stage" ]; then continue; fi
                echo "  [zip] $f"
                if have unzip; then
                    mkdir -p "$stage" && unzip -q -o "$f" -d "$stage"
                else
                    mkdir -p "$stage" && bsdtar -xf "$f" -C "$stage"
                fi
                rm -f "$f"
                changed=1
                ;;
            *.rar|*.7z)
                stage="${f%.*}"
                if [ -e "$stage" ]; then continue; fi
                echo "  [$base 后缀] $f"
                mkdir -p "$stage"
                extract_archive "$f" "$stage"
                rm -f "$f"
                changed=1
                ;;
        esac
    done < <(find "$OUTDIR" -type f -print0)

    [ "$changed" -eq 0 ] && break
    [ "$round" -gt 10 ] && { echo "警告: 递归层数超过 10, 停止" >&2; break; }
done

echo "==> 完成。文件位于 $(pwd)/$OUTDIR"
find "$OUTDIR" -type f -exec ls -lh {} + | awk '{print $5"\t"$9}'
