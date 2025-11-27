# fastp-rs

Rust 版本的 fastp 重写实现（rewrite）。本项目以高性能与可维护性为目标，复刻并增强原工具 fastp 的核心功能：

- 多线程流水线：Reader → Worker → Writer，基于通道并发
- 成对（PE）重叠纠错与统计：支持重叠区域的差异统计与平均差
- PolyX/PolyG 裁剪与质量过滤
- 报告输出：JSON 与 HTML（包含 PolyX/PolyG 与 PE Overlap 指标）
- I/O 优化：缓冲读取、批量打包、复用缓冲减少分配
- 并行压缩：内置 gzip 压缩线程池；亦支持外部 `pigz`
- 面向集群的参数调优：`pack_size`、`queue_depth`、`-w` 线程数与 `-z` 压缩等级

> 说明：本项目是对原工具 fastp 的 Rust 重写（rewrite），在不改变基本使用习惯的前提下，针对多核 CPU 和高并发 I/O 场景进行了优化。

## 上游项目链接

- 原 fastp 项目仓库：https://github.com/OpenGene/fastp

## 构建

```bash
cargo build --release
```

生成的可执行文件位于 `./target/release/fastp_rs`。

## 使用示例

- 成对（PE）输入：

```bash
./target/release/fastp_rs \
  -i R1.fq.gz -I R2.fq.gz \
  -o out1.fq.gz -O out2.fq.gz \
  --json fastp.json --html fastp.html \
  -w 24 --pack_size 30000 --queue_depth 96 -z 1 \
  -x --poly_x_min_len 10 \
  --trim_poly_g --poly_g_min_len 10 \
  -c --overlap_len_require 30 \
     --overlap_diff_limit 5 \
     --overlap_diff_percent_limit 20
```

- 单端（SE）输入：

```bash
./target/release/fastp_rs \
  -i R1.fq.gz \
  -o out1.fq.gz \
  --json fastp.json --html fastp.html \
  -w 24 --pack_size 30000 --queue_depth 96 -z 1 \
  -x --poly_x_min_len 10 \
  --trim_poly_g --poly_g_min_len 10
```

## 并行压缩

- 默认内置 gzip 压缩线程池（与处理线程协同，提高吞吐）。
- 如需使用外部 `pigz`，可启用：

```bash
./target/release/fastp_rs ... --pigz --pigz_threads 24 -o /dev/stdout -O /dev/stdout | \
tee >(pigz -p 24 > out1.fq.gz) >(pigz -p 24 > out2.fq.gz) > /dev/null
```

## Slurm 示例脚本（示例分区）

```bash
可以用 Rust 子命令生成 sbatch 脚本，减少脚本体量：

```bash
FASTP_RS_CMD=emit_sbatch FASTP_RS_PARTITION=<your_partition> FASTP_RS_CPUS=24 FASTP_RS_MEM=64G \
  ./target/release/fastp_rs --pack_size 30000 --queue_depth 96 -z 1 > run_fastp_rs.sbatch
```

生成的 `run_fastp_rs.sbatch` 包含构建与运行指令，你可以直接 `sbatch run_fastp_rs.sbatch`。
```

## 参数说明（核心）

- `-w, --thread`：工作线程数，默认取 CPU 核数
- `--pack_size`：打包大小（每批处理的记录数），增大提升吞吐但提高内存占用
- `--queue_depth`：通道队列深度，建议为 `threads * 2 ~ 4`
- `-z, --compression`：gzip 压缩等级（0~9），1 为快速；越高 CPU 开销越大
- `--pigz`、`--pigz_threads`：启用外部 pigz 并行压缩及线程数（可选）
- `--json`、`--html`：报告文件路径，HTML 包含 PolyX/PolyG 与 PE Overlap 统计
- `-x, --poly_x_min_len`、`--trim_poly_g --poly_g_min_len`：PolyX/PolyG 裁剪阈值
- `-c, --correction`、`--overlap_len_require`、`--overlap_diff_limit`、`--overlap_diff_percent_limit`：PE 重叠纠错与统计参数

## 优化与实现细节

- Reader：`BufRead::read_line` 批量读取、复用缓冲，减少 String 分配与系统调用
- Writer：实现 `std::io::Write`，支持 `write_all`；内置 gzip 压缩线程池
- 并行流水：使用 crossbeam 通道在 Reader/Workers/Writer 间传递 `Pack`
- 报告：在 HTML/JSON 中输出 PolyX/PolyG 与 PE Overlap（平均差与计数）

## 目标与兼容性

- 目标：在不改变 fastp 使用体验的前提下，针对多核与集群环境获得更高吞吐与更稳定的资源占用
- 兼容性：支持 `.fq`/`.fq.gz` 输入；PE 与 SE 模式；报告与常用裁剪/纠错参数

## 许可

本项目遵循与上游 fastp 一致的许可协议（请参见 OpenGene/fastp 仓库中的 License 声明）。
如需在本仓库中展示完整 License 内容，请在创建仓库后同步添加相同的 License 文件。
