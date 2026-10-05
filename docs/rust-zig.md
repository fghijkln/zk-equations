# Rust / Zig 曲线后端（debug^3，实验性）

v5.0 的 C 扩展回答了"纯 Python 做不到常量时间"，但发行说明里诚实地写了
"可能有 bug"。debug^3 把同样的曲线运算再用 Rust 和 Zig 各写一遍，
三个原生后端 + 纯 Python 四路差分，互相盯着。

## 为什么是 Rust 和 Zig

- **Rust**：内存安全无 GC，`u128` 原生，借用检查器在编译期就拦下了
  C 版曾经犯过的 aliasing 类 bug（这次重写时亲测：23 个 E0502）。
  是摘掉"实验性"帽子的正统路线。
- **Zig**：`zig build-lib -target aarch64-linux-android` 一行出 .so，
  不需要 NDK（自带 libc），交叉编译体验碾压。`comptime` 适合生成
  常量时间代码。

## 架构

```
core/rsext/          Rust cdylib (curve_rs)，零依赖
core/zigext/         Zig 单文件 (curve.zig)，export fn
core/nbext/ABI.md    共享 C ABI 规范（9 个 wcurve_* 函数）
core/curve_backend.py  ctypes 统一调度：rust > zig > c > python
core/curve_c.py        薄 shim，委托给 curve_backend（原有 import 不动）
core/test_backends.py  四路差分测试
```

关键设计：Rust/Zig 只暴露**纯 C ABI**（`core/nbext/ABI.md`），
Python 侧用 ctypes 调用——**不需要 Python.h**，避开了 C 扩展
`pyconfig.h` 地狱那 5 次失败构建。

后端选择：`WITNESS_CURVE_BACKEND` 环境变量，或
`curve_backend.set_backend("zig")`，默认自动选最快的可用后端。

## 验证状态（2026-10-05）

- 四路差分测试：`python3 -m core.test_backends` → ALL BACKENDS AGREE
- 全量测试：102/102 通过（默认 rust 后端）
- Zig Android .so：本地交叉编译验证，ARM64 ELF 有效
- Rust Android：`cargo check --target aarch64-linux-android` 通过，
  最终链接在 CI 用 NDK clang 完成

## 诚实声明

- 三个原生实现都没审计。差分测试只能证明"和 Python 版一致"，
  不能证明"没有 bug"——如果 Python 版有 bug，四个一起错。
- 常量时间：操作序列是固定的，但编译器可能引入分支，没看汇编。
- 这是 debug 版，稳定版 v5.0 不受影响。
