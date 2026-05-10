# 🔧 Quarto 网站修复总结

## 问题诊断

你的网站有两个主要问题：

### 1. ❌ 资源路径问题
**问题**: 图片、视频和数据文件无法显示  
**原因**: `docs/` 子目录中的 `.qmd` 文件使用了相对于项目根目录的路径（如 `Figure/`），但生成的 HTML 在 `_build/html/docs/` 中，所以需要使用 `../Figure/` 来指向上一级目录。

**已修复**: ✅ 所有 10 个章节文件中的路径已从：
- `src="Figure/` → `src="../Figure/`
- `src="Data/` → `src="../Data/`

### 2. ❌ 资源复制配置
**问题**: `Figure/` 和 `Data/` 文件夹没有被明确配置为要复制的资源

**已修复**: ✅ 在 `_quarto.yml` 中添加了 `resources:` 配置

---

## 已应用的修复

### 修复 1: 文件路径 (10 个文件)
```
docs/chap0_preface.qmd
docs/chap1.qmd
docs/chap2.qmd
docs/chap3.qmd
docs/chap4.qmd
docs/chap5.qmd
docs/chap6.qmd
docs/chap7.qmd
docs/chap8.qmd
docs/chap10_appendix.qmd
```

**示例修改**:
```markdown
# 之前
<img src="Figure/figure_fft_time_frequency.png" />

# 之后
<img src="../Figure/figure_fft_time_frequency.png" />
```

### 修复 2: _quarto.yml 配置
```yaml
project:
  resources:
    - Figure/
    - Data/
```

---

## 现在你需要做什么

### 步骤 1: 清理旧的构建文件
```bash
cd "D:\OneDrive\Documents\GitHub\Time-Series-Data-Analysis"
rm -rf _build
```

### 步骤 2: 重新构建网站
```bash
quarto render
```

### 步骤 3: 预览修复后的网站
```bash
quarto preview
```

然后在浏览器中访问 http://localhost:4200

---

## 验证修复

检查以下内容是否都能正常显示：

✅ **首页** - `index.qmd`
- [ ] 所有示例图片显示正确
- [ ] 导航菜单工作正常

✅ **第 0 章 (Preface)** - `chap0_preface.qmd`
- [ ] 所有图片显示
- [ ] 音频播放器工作
- [ ] 标签页切换正常

✅ **其他章节**
- [ ] 图表和图片显示
- [ ] 代码块显示

---

## 如果仍有问题

### 问题: 图片仍然显示不出来

**解决方案 1**: 检查浏览器控制台（F12）
- 打开开发者工具
- 查看 Console 标签页
- 检查是否有错误信息

**解决方案 2**: 检查文件是否存在
```bash
ls -la Figure/
ls -la Data/
```

### 问题: 某些链接还是打不开

**可能的原因**:
1. 内部链接使用了 docsify 格式（如 `[link](chap1.md)`）
2. 需要改为 Quarto 格式（如 `[link](chap1.qmd)` 或 `[link](../docs/chap1.qmd)`）

**修复方法**: 
- 使用 `.qmd` 替代 `.md` 扩展名
- 使用相对路径指向其他文档

---

## 完整清单

- ✅ 修复了所有资源路径（Figure/ 和 Data/）
- ✅ 更新了 _quarto.yml 资源配置
- ✅ 创建了修复脚本 (`fix_paths.py`)
- 🔄 **下一步**: 在本地重新构建网站

---

## 后续改进建议

1. **使用 Markdown 图片语法** 替代 HTML：
   ```markdown
   # 好的做法
   ![alt text](../Figure/image.png)
   
   # 避免
   <img src="..." />
   ```

2. **使用相对链接** 链接到其他页面：
   ```markdown
   # 好的做法
   [链接到第1章](chap1.qmd)
   
   # 避免
   [链接到第1章](../docs/chap1.qmd)
   ```

3. **定期检查** 生成的 HTML 是否有断掉的图片

---

**修复完成时间**: 2026-05-09

需要帮助？查看 `QUICKSTART.md` 或 `MIGRATION.md`！
