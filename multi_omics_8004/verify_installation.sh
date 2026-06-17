#!/bin/bash
# Verification script for Multi-Omics Analysis Platform

echo "=========================================="
echo "Multi-Omics Platform Installation Check"
echo "=========================================="
echo ""

# Check Python version
echo "1. Checking Python version..."
python3 --version
echo ""

# Check if all required files exist
echo "2. Checking required files..."
files=(
    "app.py"
    "requirements.txt"
    "config/llm_config.json"
    "src/executor.py"
    "src/planner.py"
    "src/validator.py"
    "src/xai_explainer.py"
    "src/arbiter.py"
    "src/registry.py"
    "src/models.py"
    "src/llm_client.py"
    "src/rag_retriever.py"
    "templates/index.html"
)

all_exist=true
for file in "${files[@]}"; do
    if [ -f "$file" ]; then
        echo "  ✓ $file"
    else
        echo "  ✗ $file (MISSING)"
        all_exist=false
    fi
done
echo ""

# Check if directories exist
echo "3. Checking required directories..."
dirs=("data" "output" "logs" "static" "templates" "src" "config")
for dir in "${dirs[@]}"; do
    if [ -d "$dir" ]; then
        echo "  ✓ $dir/"
    else
        echo "  ✗ $dir/ (MISSING)"
        all_exist=false
    fi
done
echo ""

# Check API key in config
echo "4. Checking API key..."
if grep -q "YOUR_API_KEY_HERE" config/llm_config.json; then
    echo "  ✓ API key updated correctly"
else
    echo "  ✗ API key not updated"
    all_exist=false
fi
echo ""

# Check if app.py has real analysis functions
echo "5. Checking app.py for real analysis functions..."
if grep -q "compute_volcano_data" app.py && grep -q "compute_umap_data" app.py; then
    echo "  ✓ Real analysis functions present"
else
    echo "  ✗ Real analysis functions missing"
    all_exist=false
fi
echo ""

# Summary
echo "=========================================="
if [ "$all_exist" = true ]; then
    echo "✓ All checks passed!"
    echo "Ready for deployment."
else
    echo "✗ Some checks failed."
    echo "Please review the errors above."
fi
echo "=========================================="
