# 🚀 构建指南 - 修复所有问题

## ✅ 已完成的修复

1. ✅ **恢复英文内容** - 所有页面现在都是英文
2. ✅ **修复图片路径** - 从 `../Figure/` 改为 `Figure/`
3. ✅ **修复音频路径** - 从 `../Data/` 改为 `Data/`
4. ✅ **更新资源配置** - `_quarto.yml` 现在明确指定资源目录

## 🏗️ 构建步骤

### 在 WSL 中执行：

```bash
cd /mnt/d/OneDrive/Documents/GitHub/Time-Series-Data-Analysis

# 1. 完全清空旧构建
rm -rf _build

# 2. 清理缓存
rm -rf _freeze

# 3. 全新构建（这会复制所有资源）
quarto render

# 4. 启动预览
quarto preview
```

### 或者一条命令：

```bash
cd /mnt/d/OneDrive/Documents/GitHub/Time-Series-Data-Analysis && \
rm -rf _build _freeze && \
quarto render && \
quarto preview
```

## ✅ 验证资源是否正确复制

构建完成后，在 WSL 中检查：

```bash
# 检查 Figure 文件夹是否被复制到输出目录
ls _build/html/Figure/ | head -5

# 检查 Data 文件夹是否被复制到输出目录
ls _build/html/Data/ | head -5
```

**预期输出**：
- `Figure/` 目录中应该有所有图片文件
- `Data/` 目录中应该有所有音频/视频文件

## 📋 路径配置说明

**在 `.qmd` 文件中引用资源的正确方式**：

```markdown
# 图片
![Description](Figure/image_name.png)

# 或使用 HTML
<img src="Figure/image_name.png" width="100%"/>

# 音频
<audio controls>
  <source src="Data/audio_file.wav" type="audio/wav">
</audio>

# 视频
<video controls width="100%">
  <source src="Data/video_file.mp4" type="video/mp4">
</video>
```

## 🔍 如果图片/音频仍然无法显示

### 检查 1：验证文件存在

```bash
# 检查源文件是否存在
ls -la Figure/ | grep -E "\.png|\.jpg"
ls -la Data/ | grep -E "\.wav|\.mp4|\.mov"
```

### 检查 2：验证输出目录结构

```bash
# 检查构建输出中是否有资源
ls -la _build/html/

# 应该看到：
# Figure/
# Data/
# index.html
# chap0_preface.html
# 等等
```

### 检查 3：在浏览器中验证

1. 打开 http://localhost:4200
2. 按 F12 打开开发者工具
3. 查看 Console 标签页，看是否有错误信息
4. 查看 Network 标签页，检查资源请求是否返回 404

### 检查 4：查看源代码

在浏览器中右键 → "查看页面源代码"，寻找类似的链接：
```html
<img src="Figure/figure_fft_time_frequency.png" ...>
<audio src="Data/indian_cuckoo_short.wav" ...>
```

## 🎯 最终测试清单

- [ ] `quarto render` 完成无错误
- [ ] `_build/html/Figure/` 目录存在且包含图片
- [ ] `_build/html/Data/` 目录存在且包含音频/视频
- [ ] http://localhost:4200 能打开
- [ ] 首页（index）的图片和音频都能显示
- [ ] chap0_preface 的所有图片和音频都能播放
- [ ] 其他章节的内容正常显示

## 💡 提示

如果还有问题，可以尝试：

```bash
# 清理所有 Quarto 缓存
quarto clean

# 再次构建
quarto render

# 使用详细输出查看发生了什么
quarto render --verbose
```

## 📝 现在就开始吧！

```bash
cd /mnt/d/OneDrive/Documents/GitHub/Time-Series-Data-Analysis && \
rm -rf _build _freeze && \
quarto render && \
quarto preview
```

然后在浏览器中打开 http://localhost:4200，享受你的网站！ 🎉
