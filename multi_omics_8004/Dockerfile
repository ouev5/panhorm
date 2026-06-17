# 多阶段构建：builder阶段
FROM python:3.13-slim AS builder

# 设置环境变量
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

# 安装系统依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
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
    && rm -rf /var/lib/apt/lists/*

# 创建并设置工作目录
WORKDIR /app

# 复制依赖文件
COPY requirements.txt requirements-dev.txt ./

# 安装Python依赖
RUN pip install --upgrade pip \
    && pip install -r requirements.txt \
    && pip install -r requirements-dev.txt

# 复制项目文件
COPY . .

# 多阶段构建：runtime阶段
FROM python:3.13-slim AS runtime

# 设置环境变量
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    APP_ENV=production \
    DATA_DIR=/app/data \
    OUTPUT_DIR=/app/output \
    LOGS_DIR=/app/logs

# 安装系统依赖（运行时需要的）
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    git \
    libz-dev \
    && rm -rf /var/lib/apt/lists/*

# 创建非root用户
RUN groupadd -g 1000 appuser \
    && useradd -u 1000 -g appuser -s /bin/bash -m appuser

# 创建必要的目录
RUN mkdir -p $DATA_DIR $OUTPUT_DIR $LOGS_DIR \
    && chown -R appuser:appuser /app

# 切换到非root用户
USER appuser

# 设置工作目录
WORKDIR /app

# 从builder阶段复制依赖
COPY --from=builder /usr/local/lib/python3.13/site-packages /usr/local/lib/python3.13/site-packages
COPY --from=builder /usr/local/bin /usr/local/bin

# 复制项目文件
COPY --chown=appuser:appuser . .

# 健康检查端点
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD curl -f http://localhost:7860/ || exit 1

# 暴露端口
EXPOSE 7860

# 启动命令
CMD ["python", "src/webui/app.py"]
