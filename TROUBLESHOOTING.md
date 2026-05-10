# 🔍 故障排除指南

## 问题：重命名错误 (rename error)

```
ERROR: NotFound: No such file or directory (os error 2): 
rename '/path/to/docs/chap0_preface.html' -> '/path/to/_build/html/docs/chap0_preface.html'
```

### 原因
Quarto 在尝试将渲染的文件移动到输出目录时，目标目录不存在。

### 解决方案

#### 方法 1: 手动清理并重新构建 (推荐)

```bash
cd "D:\OneDrive\Documents\GitHub\Time-Series-Data-Analysis"

# 1. 完全清理构建文件夹
rmdir /s /q _build
# 或 macOS/Linux:
# rm -rf _build

# 2. 重新构建
quarto render

# 3. 预览
quarto preview
```

#### 方法 2: 使用 Quarto 的 clean 选项

```bash
cd "D:\OneDrive\Documents\GitHub\Time-Series-Data-Analysis"

# 清理并重新渲染
quarto render --execute
```

#### 方法 3: 检查磁盘空间和权限

如果上述方法都不行，可能是权限问题：

```bash
# 检查 _build 目录权限
ls -la _build/

# 检查磁盘空间
df -h

# 检查文件是否被其他程序占用
# Windows: 检查资源管理器中是否打开了 _build 目录
# macOS/Linux: lsof | grep _build
```

---

## 常见问题

### Q1: 网站渲染后图片仍然显示不出来

**检查清单**:
1. 打开浏览器开发者工具 (F12)
2. 查看 Network 标签页
3. 找到失败的图片请求
4. 检查 URL 是否正确

**常见原因和修复**:
```
问题: 图片 URL 为 /Figure/image.png
修复: 应该是 ../Figure/image.png (如果在 docs 子目录中)

问题: 文件大小写不匹配 (Windows 忽略大小写，Linux 不忽略)
修复: 确保文件名大小写与代码中完全相同
```

### Q2: 某些链接打不开

**检查**:
- [ ] 链接是否使用了正确的 `.qmd` 扩展名
- [ ] 链接的相对路径是否正确
- [ ] 文件是否真的存在

**示例**:
```markdown
# ✅ 正确
[链接](chap1.qmd)
[链接](../docs/chap2.qmd)

# ❌ 错误 (Docsify 格式，Quarto 不支持)
[链接](chap1.md)
[链接](#/docs/chap1)
```

### Q3: 网站样式看起来不对

**可能原因**:
1. `styles.css` 没有被加载
2. Quarto 主题配置问题

**修复**:
```yaml
# 在 _quarto.yml 中确保有正确的主题配置
format:
  html:
    theme:
      light: cosmo
      dark: darkly
    css: styles.css
```

---

## 深度故障排除

### 如果重新构建仍然失败

#### 步骤 1: 检查 Quarto 安装

```bash
quarto --version
quarto check
```

输出应该显示：
```
Quarto 1.3.x or newer
[✓] Checking Quarto installation...
```

#### 步骤 2: 验证文件结构

```bash
# 确保所需的文件都存在
ls -la index.qmd
ls -la _quarto.yml
ls -la styles.css
ls -la docs/chap*.qmd
ls -la Figure/
ls -la Data/
```

#### 步骤 3: 尝试渲染单个文件

```bash
# 渲染主页看看是否有效
quarto render index.qmd

# 如果主页可以，尝试渲染一个章节
quarto render docs/chap0_preface.qmd
```

#### 步骤 4: 查看详细的错误信息

```bash
# 启用详细日志
quarto render --verbose
```

---

## 最终解决方案

如果以上方法都不行，尝试这个**完整重置**:

```bash
cd "D:\OneDrive\Documents\GitHub\Time-Series-Data-Analysis"

# 1. 删除所有构建文件
rmdir /s /q _build          # Windows
# rm -rf _build             # macOS/Linux

# 2. 删除 docs 目录中的临时文件
cd docs
del *.html                  # Windows
# rm -f *.html              # macOS/Linux
cd ..

# 3. 重新渲染整个项目
quarto render

# 4. 启动预览服务器
quarto preview
```

---

## 如果问题仍未解决

请提供以下信息进行诊断：

1. Quarto 版本: `quarto --version`
2. 操作系统版本
3. 完整的错误消息
4. `_quarto.yml` 的内容
5. 项目文件夹的结构: `tree` 或 `dir /s`

---

## 相关资源

- [Quarto 官方文档](https://quarto.org/docs/)
- [Quarto 网站项目](https://quarto.org/docs/projects/quarto-projects.html)
- [Quarto 故障排除](https://quarto.org/docs/getting-started/authoring.html)

---

**最后更新**: 2026-05-09

有问题？查看 `QUICKSTART.md` 或 `FIX_SUMMARY.md`！
