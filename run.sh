#!/bin/bash

# 简化版执行脚本
set -e  # 遇到错误立即退出

echo "=== 开始执行项目编译和实验 ==="

# 步骤1: 在 wakson 目录编译
echo "1. 在 /home/project/wakson/ 下编译..."
cd /home/project/wakson

echo "执行 make clean..."
make clean || echo "警告: make clean 可能有误，但继续..."

echo "执行 make..."
make 

echo "✓ 编译完成"

# 步骤2: 运行实验脚本
echo -e "\n2. 在 Application 目录运行实验..."
cd /home/project/wakson/Application

echo "检查并运行 run_experiments.py..."

if [ ! -f "run_experiments.py" ]; then
    echo "错误: 未找到 run_experiments.py"
    exit 1
fi

# 给脚本添加执行权限（如果需要）
chmod +x run_experiments.py 2>/dev/null || true

# 运行脚本
if [[ $(head -n1 run_experiments.py) == "#!"* ]]; then
    ./run_experiments.py
else
    python3 run_experiments.py
fi

echo "✓ 所有步骤完成！"