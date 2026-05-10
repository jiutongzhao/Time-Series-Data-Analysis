# Docsify → Quarto 迁移指南

## 完成的迁移任务

✅ 创建了 `_quarto.yml` 配置文件
✅ 将所有 markdown 文件转换为 `.qmd` 格式
✅ 创建了 `index.qmd` 主页
✅ 创建了自定义 `styles.css` 样式
✅ 配置了 KaTeX 数学支持
✅ 配置了 Mermaid 图表支持

## 文件结构

```
Time-Series-Data-Analysis/
├── _quarto.yml              # Quarto 项目配置
├── index.qmd               # 主页
├── styles.css              # 自定义样式
├── docs/
│   ├── chap0_preface.qmd
│   ├── chap1.qmd
│   ├── chap2.qmd
│   ├── ... 其他章节
│   └── chap10_appendix.qmd
├── Figure/                 # 图片文件夹（保持不变）
├── Data/                   # 数据文件夹（保持不变）
└── _build/                 # 构建输出文件夹（自动生成）
    └── html/               # HTML 网站输出
```

## 如何构建网站

### 前置要求

1. **安装 Quarto** (https://quarto.org)
   - Windows: 下载安装程序或使用 `choco install quarto`
   - macOS: 使用 `brew install quarto`
   - Linux: 参考官方文档

2. **验证安装**
   ```bash
   quarto --version
   ```

### 构建步骤

1. **进入项目目录**
   ```bash
   cd D:\OneDrive\Documents\GitHub\Time-Series-Data-Analysis
   ```

2. **渲染网站**
   ```bash
   quarto render
   ```

3. **预览网站**
   ```bash
   quarto preview
   ```
   这会在浏览器中自动打开 `http://localhost:4200`

### GitHub Pages 部署

如果你想在 GitHub Pages 上部署：

1. **创建 GitHub Actions 工作流**

   创建文件 `.github/workflows/publish.yml`:
   ```yaml
   name: Publish Website
   
   on:
     push:
       branches:
         - main
     pull_request: {}
   
   permissions:
     contents: read
     pages: write
     id-token: write
   
   concurrency:
     group: "pages"
     cancel-in-progress: false
   
   jobs:
     build-deploy:
       runs-on: ubuntu-latest
       steps:
         - uses: actions/checkout@v4
         
         - name: Set up Quarto
           uses: quarto-dev/quarto-actions/setup@v2
         
         - name: Render Quarto project
           uses: quarto-dev/quarto-actions/render@v2
         
         - name: Upload pages artifact
           uses: actions/upload-pages-artifact@v3
           with:
             path: '_build/html'
         
         - name: Deploy to GitHub Pages
           id: deployment
           uses: actions/deploy-pages@v4
   ```

2. **在 GitHub 设置中启用 Pages**
   - 进入仓库 Settings → Pages
   - Source 选择 "GitHub Actions"

## 迁移更改说明

### Docsify → Quarto 的转换

| 功能 | Docsify | Quarto |
|------|---------|--------|
| 配置文件 | `index.html` | `_quarto.yml` |
| 导航 | `_sidebar.md` | `_quarto.yml` 中的 sidebar |
| 标签页 | `<!-- tabs:start -->` | `::: {.panel-tabset}` |
| 输出目录 | `docs/` | `_build/html/` |
| 主题 | CDN 主题 | Quarto 内置主题 |

### 保持不变

- ✅ 所有原始 markdown 内容
- ✅ 所有图片和媒体文件位置
- ✅ 数学公式（从 KaTeX 到 Quarto 的 KaTeX）
- ✅ Mermaid 图表

## 主要优势

1. **更好的 R/Python 集成**: 支持在 Quarto 文档中执行代码块
2. **多格式输出**: 轻松生成 HTML、PDF、Word 等多种格式
3. **改进的数学支持**: 原生 LaTeX 和 MathJax 支持
4. **现代工具链**: 基于 Pandoc 的可靠转换
5. **社区支持**: 活跃的开发和广泛的文档

## 常见问题

### Q: 如何修改网站外观?
A: 编辑 `_quarto.yml` 中的 `format.html.theme` 部分。可用主题有: `cosmo`, `darkly`, `journal`, `lumen`, `quartz` 等。

### Q: 如何添加新页面?
A: 
1. 在 `docs/` 中创建新的 `.qmd` 文件
2. 在 `_quarto.yml` 的 `sidebar.contents` 中添加引用

### Q: 原始 docsify 文件现在用不了吗?
A: 原始 `.md` 文件和 `index.html` 仍然保留，但不再使用。如果需要保留，可以保留在版本控制中。

### Q: 如何获得更好的性能?
A: Quarto 已经优化了构建速度。使用 `quarto render --no-execute` 可以跳过代码执行（如果有 Python 代码块）。

## 故障排除

### 网站不显示图片?
- 检查图片路径是否相对于项目根目录
- 确保 `Figure/` 文件夹存在且内容完整

### 数学公式显示有问题?
- 检查 `_quarto.yml` 中的 `html-math-method: katex` 配置
- 确保使用正确的 LaTeX 语法

### 构建失败?
- 运行 `quarto check` 诊断环境问题
- 查看错误消息中的具体文件和行号
- 确保所有 markdown 文件都有有效的 YAML frontmatter（或不需要）

## 后续步骤

1. 本地测试网站: `quarto preview`
2. 提交更改到 GitHub
3. 设置 GitHub Pages 自动部署（见上文）
4. 更新项目 README，指向新的 Quarto 网站

---

**迁移完成日期**: 2026-05-09
