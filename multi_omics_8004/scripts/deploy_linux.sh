#!/bin/bash

# 部署脚本 (Linux)

set -e

echo "开始部署 autoBA_optimized..."

# 定义变量
DEPLOY_DIR="$(pwd)"
PROJECT_DIR="/opt/autoBA_optimized"
VENV_DIR="$PROJECT_DIR/venv"
SERVICE_NAME="autoBA.service"
NGINX_CONF="/etc/nginx/sites-available/autoBA"
NGINX_ENABLED="/etc/nginx/sites-enabled/autoBA"

# 1. 安装依赖
echo "安装系统依赖..."
sudo apt update
sudo apt install -y \
    build-essential \
    curl \
    wget \
    git \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libffi-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    zlib1g-dev \
    nginx \
    supervisor

# 2. 安装 Python 3.13
echo "安装 Python 3.13..."

# 检查是否已安装 Python 3.13
if ! command -v python3.13 &> /dev/null; then
    # 下载并安装 Python 3.13
    wget https://www.python.org/ftp/python/3.13.0/Python-3.13.0.tar.xz
    tar -xf Python-3.13.0.tar.xz
    cd Python-3.13.0
    ./configure --enable-optimizations
    make -j $(nproc)
    sudo make altinstall
    cd ..
    rm -rf Python-3.13.0 Python-3.13.0.tar.xz
    echo "Python 3.13 安装完成"
else
    echo "Python 3.13 已安装"
fi

# 3. 创建项目目录
echo "创建项目目录..."
sudo mkdir -p $PROJECT_DIR
sudo chown $USER:$USER $PROJECT_DIR

# 4. 解压项目文件
echo "解压项目文件..."
if [ -f "autoBA_optimized.zip" ]; then
    unzip autoBA_optimized.zip -d $PROJECT_DIR
else
    echo "错误：autoBA_optimized.zip 文件不存在"
    exit 1
fi

# 5. 创建虚拟环境
echo "创建虚拟环境..."
python3.13 -m venv $VENV_DIR

# 6. 安装Python依赖
echo "安装Python依赖..."
$VENV_DIR/bin/pip install --upgrade pip
if [ -f "requirements-linux.txt" ]; then
    $VENV_DIR/bin/pip install -r requirements-linux.txt
else
    $VENV_DIR/bin/pip install -r $PROJECT_DIR/requirements.txt
fi

# 7. 创建systemd服务
echo "创建systemd服务..."
SERVICE_CONTENT="[Unit]
Description=autoBA_optimized Service
After=network.target

[Service]
User=$USER
WorkingDirectory=$PROJECT_DIR
ExecStart=$VENV_DIR/bin/python src/webui/app.py
Restart=on-failure
RestartSec=5s
Environment=APP_ENV=production
Environment=LOG_LEVEL=INFO
Environment=MAX_CONCURRENT_TASKS=4
Environment=DATA_DIR=$PROJECT_DIR/data
Environment=OUTPUT_DIR=$PROJECT_DIR/output
Environment=LOGS_DIR=$PROJECT_DIR/logs

[Install]
WantedBy=multi-user.target"

sudo bash -c "echo '$SERVICE_CONTENT' > /etc/systemd/system/$SERVICE_NAME"

# 8. 创建必要的目录
echo "创建必要的目录..."
mkdir -p $PROJECT_DIR/data $PROJECT_DIR/output $PROJECT_DIR/logs $PROJECT_DIR/config/pipelines

# 9. 配置Nginx反向代理
echo "配置Nginx反向代理..."
NGINX_CONTENT="server {
    listen 80;
    server_name _;

    location / {
        proxy_pass http://localhost:7860;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }

    error_log /var/log/nginx/autoBA_error.log;
    access_log /var/log/nginx/autoBA_access.log;
}"

sudo bash -c "echo '$NGINX_CONTENT' > $NGINX_CONF"

# 启用Nginx配置
sudo ln -sf $NGINX_CONF $NGINX_ENABLED

# 10. 启动服务
echo "启动服务..."

# 重新加载systemd
sudo systemctl daemon-reload

# 启动autoBA服务
sudo systemctl start $SERVICE_NAME

# 启用服务自启动
sudo systemctl enable $SERVICE_NAME

# 重启Nginx
sudo systemctl restart nginx

# 11. 验证部署
echo "验证部署..."
sleep 5

# 检查服务状态
echo "服务状态："
sudo systemctl status $SERVICE_NAME

# 检查Nginx状态
echo "Nginx状态："
sudo systemctl status nginx

# 检查端口
echo "检查端口："
netstat -tuln | grep 7860

# 完成
echo "\n部署完成！"
echo "Web界面地址：http://$(hostname -I | awk '{print $1}')"
echo "服务状态可通过以下命令查看：sudo systemctl status $SERVICE_NAME"
echo "日志可通过以下命令查看：journalctl -u $SERVICE_NAME"
