# 🚀 快速开始 - Quarto 网站

## 第一步：安装 Quarto

如果你还没安装 Quarto，请先安装：

**Windows:**
```bash
# 使用 Chocolatey
choco install quarto

# 或者下载安装程序
# https://quarto.org/docs/get-started/
```

**macOS:**
```bash
brew install quarto
```

**Linux:**
参考 https://quarto.org/docs/get-started/

验证安装：
```bash
quarto --version
```

## 第二步：生成网站

在项目根目录执行：

```bash
cd "D:\OneDrive\Documents\GitHub\Time-Series-Data-Analysis"

# 渲染网站（生成 HTML）
quarto render

# 或者直接预览（自动生成并在浏览器打开）
quarto preview
```

## 第三步：访问网站

- **本地预览**: `quarto preview` 后会自动打开 http://localhost:4200
- **直接打开**: 在浏览器打开 `_build/html/index.html`

## 常见问题排查

### 问题 1: `quarto: command not found`

**解决方案**: 
- 检查 Quarto 是否已安装: `quarto --version`
- 如果未安装，请参考第一步
- 如果已安装但找不到，重启终端或重启电脑

### 问题 2: 找不到图片或资源文件

**解决方案**:
- 确保 `Figure/` 和 `Data/` 文件夹在项目根目录
- 检查文件路径是否正确（相对于项目根目录）
- 如果在本地预览，确保文件访问权限正确

### 问题 3: 数学公式显示有问题

**解决方案**:
- 检查 LaTeX 语法是否正确
- 确保 `_quarto.yml` 中有 `html-math-method: katex`
- 刷新浏览器缓存 (Ctrl+Shift+Delete)

## 文件结构

```
Time-Series-Data-Analysis/
├── index.qmd                    # 主页
├── _quarto.yml                  # Quarto 配置 ⭐ 重要
├── styles.css                   # 自定义样式
├── docs/
│   ├── chap0_preface.qmd        # 序言
│   ├── chap1.qmd                # 第 1 章
│   ├── ... 其他章节
│   └── chap10_appendix.qmd      # 附录
├── Figure/                      # 图片文件夹
├── Data/                        # 数据文件夹
└── _build/                      # 构建输出（自动生成）
    └── html/                    # 最终的 HTML 网站
```

## 部署到 GitHub Pages

如果你想让网站在线访问，可以部署到 GitHub Pages：

### 方式 1: 自动部署（推荐）

在 `.github/workflows/publish.yml` 中配置自动构建和部署。参考 `MIGRATION.md` 中的详细步骤。

### 方式 2: 手动部署

1. 本地生成网站:
   ```bash
   quarto render
   ```

2. 上传 `_build/html/` 文件夹到 GitHub 的 `gh-pages` 分支：
   ```bash
   git add _build/
   git commit -m "Build website"
   git push origin main
   ```

3. 在 GitHub 仓库设置中启用 Pages，指向 `gh-pages` 分支

## 编辑文档

要添加或修改内容：

1. **编辑现有文件**: 直接修改 `.qmd` 文件
2. **添加新章节**: 
   - 在 `docs/` 中创建新的 `.qmd` 文件
   - 在 `_quarto.yml` 的 sidebar 中添加引用

例如，添加新章节 `chap9.qmd`:
```yaml
# 在 _quarto.yml 中
sidebar:
  - id: docs
    contents:
      - section: "New Section"
        contents:
          - docs/chap9.qmd
```

## 有用的命令

```bash
# 检查 Quarto 环境
quarto check

# 生成网站（不自动打开）
quarto render

# 实时预览（自动刷新）
quarto preview

# 生成特定格式
quarto render --to html
quarto render --to pdf      # 需要 TinyTeX

# 清理构建文件
rm -rf _build
```

## 更多信息

- 📖 [Quarto 官方文档](https://quarto.org)
- 🎨 [可用主题列表](https://quarto.org/docs/output-formats/html-themes.html)
- 📝 [Markdown 语法](https://quarto.org/docs/authoring/markdown-basics.html)
- 🔧 [高级配置](https://quarto.org/docs/projects/quarto-projects.html)

---

需要帮助？查看 `MIGRATION.md` 获取更多详细信息！
