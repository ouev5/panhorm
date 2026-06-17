# 部署准备脚本 (Windows)

# 设置错误时停止执行
$ErrorActionPreference = "Stop"

# 定义变量
$projectRoot = Split-Path -Parent $PSScriptRoot
$outputDir = Join-Path $projectRoot "deploy"
$dockerImageName = "autoba/optimized:latest"

# 创建输出目录
if (-not (Test-Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir -Force
    Write-Host "创建输出目录: $outputDir"
}

# 清理旧文件
Write-Host "清理旧文件..."
Remove-Item -Path "$outputDir\*" -Recurse -Force -ErrorAction SilentlyContinue

# 1. 打包项目文件
Write-Host "打包项目文件..."
$zipFile = Join-Path $outputDir "autoBA_optimized.zip"
Compress-Archive -Path "$projectRoot\*" -DestinationPath $zipFile -Force
Write-Host "项目文件已打包到: $zipFile"

# 2. 生成Linux依赖文件
Write-Host "生成Linux依赖文件..."
$requirementsLinux = Join-Path $outputDir "requirements-linux.txt"

# 读取原始依赖文件
$requirements = Get-Content (Join-Path $projectRoot "requirements.txt")

# 处理依赖，确保兼容Linux
$linuxRequirements = @()
foreach ($line in $requirements) {
    # 跳过注释和空行
    if ($line -match '^\s*#|^\s*$') {
        continue
    }
    # 添加到Linux依赖文件
    $linuxRequirements += $line
}

# 写入Linux依赖文件
$linuxRequirements | Out-File -FilePath $requirementsLinux -Encoding UTF8
Write-Host "Linux依赖文件已生成: $requirementsLinux"

# 3. 构建Docker镜像（跨平台）
Write-Host "构建Docker镜像..."
try {
    # 检查Docker是否安装
    docker --version
    
    # 切换到项目根目录
    Push-Location $projectRoot
    
    # 构建镜像
    docker build --platform linux/amd64 -t $dockerImageName .
    
    # 保存镜像到文件
    $dockerImageFile = Join-Path $outputDir "autoBA_optimized.tar"
    docker save -o $dockerImageFile $dockerImageName
    
    Write-Host "Docker镜像已构建并保存到: $dockerImageFile"
    
    # 清理Docker镜像
    docker rmi $dockerImageName
    
    Pop-Location
} catch {
    Write-Host "构建Docker镜像失败: $($_.Exception.Message)"
    Write-Host "请确保Docker已安装并运行"
}

# 4. 复制部署脚本
Write-Host "复制部署脚本..."
$deployScript = Join-Path $outputDir "deploy_linux.sh"
Copy-Item -Path (Join-Path $PSScriptRoot "deploy_linux.sh") -Destination $deployScript

# 设置脚本执行权限（在Linux端执行）
Write-Host "设置部署脚本执行权限..."
# 注意：在Windows上无法直接设置Linux权限，需要在Linux端执行chmod

# 5. 生成部署文档
Write-Host "生成部署文档..."
$deployDoc = Join-Path $outputDir "DEPLOYMENT.md"
$docContent = @'
# autoBA_optimized 部署指南

## 1. 准备工作

### 系统要求
- Linux服务器（推荐Ubuntu 22.04+）
- Docker 20.10+
- Python 3.13
- Nginx
- 至少4GB内存
- 至少20GB磁盘空间

### 部署包内容
- `autoBA_optimized.zip` - 项目源代码
- `autoBA_optimized.tar` - Docker镜像
- `requirements-linux.txt` - Linux依赖文件
- `deploy_linux.sh` - Linux部署脚本

## 2. 部署步骤

### 方法一：使用Docker部署（推荐）

1. **上传部署包到服务器**
   ```bash
   scp autoBA_optimized.tar user@server:/path/to/deploy/
   ```

2. **加载Docker镜像**
   ```bash
   docker load -i autoBA_optimized.tar
   ```

3. **运行容器**
   ```bash
   docker run -d --name autoBA \
     -p 7860:7860 \
     -v ./data:/app/data \
     -v ./output:/app/output \
     -v ./logs:/app/logs \
     -v ./config:/app/config \
     autoba/optimized:latest
   ```

### 方法二：使用系统服务部署

1. **上传部署包到服务器**
   ```bash
   scp autoBA_optimized.zip user@server:/path/to/deploy/
   scp deploy_linux.sh user@server:/path/to/deploy/
   scp requirements-linux.txt user@server:/path/to/deploy/
   ```

2. **执行部署脚本**
   ```bash
   chmod +x deploy_linux.sh
   ./deploy_linux.sh
   ```

## 3. 验证部署

1. **访问Web界面**
   打开浏览器，访问 `http://your-server-ip:7860`

2. **检查服务状态**
   - Docker部署：`docker ps`
   - 系统服务：`systemctl status autoBA.service`

3. **检查日志**
   - Docker部署：`docker logs autoBA`
   - 系统服务：`journalctl -u autoBA.service`

## 4. 配置

### 环境变量
- `APP_ENV` - 应用环境（production/development）
- `LOG_LEVEL` - 日志级别（INFO/WARNING/ERROR）
- `MAX_CONCURRENT_TASKS` - 最大并发任务数
- `DATA_DIR` - 数据目录
- `OUTPUT_DIR` - 输出目录
- `LOGS_DIR` - 日志目录

### Nginx配置
（如果使用方法二部署）

## 5. 故障排除

### 常见问题
1. **服务无法启动**
   - 检查端口是否被占用
   - 检查Docker是否运行
   - 检查系统服务状态

2. **Web界面无法访问**
   - 检查防火墙设置
   - 检查Nginx配置
   - 检查服务日志

3. **分析任务失败**
   - 检查输入数据格式
   - 检查工具依赖是否安装
   - 检查资源使用情况

## 6. 维护

### 升级
1. **停止服务**
   - Docker：`docker stop autoBA && docker rm autoBA`
   - 系统服务：`systemctl stop autoBA.service`

2. **更新部署包**
   - 上传新的部署包
   - 重复部署步骤

### 备份
1. **备份数据**
   ```bash
   tar -czf autoBA_backup_$(date +%Y%m%d).tar.gz data/ output/ config/
   ```

2. **恢复备份**
   ```bash
   tar -xzf autoBA_backup_*.tar.gz
   ```
'@

$docContent | Out-File -FilePath $deployDoc -Encoding UTF8
Write-Host "部署文档已生成: $deployDoc"

# 完成
Write-Host "\n部署准备完成！"
Write-Host "部署包已生成到: $outputDir"
Write-Host "请参考部署文档进行部署。"
