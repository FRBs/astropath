#!/bin/bash
# Verification script for assign_host documentation

echo "=========================================="
echo "Documentation Verification"
echo "=========================================="

# Check RST files exist
echo -e "\n1. Checking RST files..."
for file in assign_host.rst quickstart_assign_host.rst; do
    if [ -f "_build/html/${file%.rst}.html" ]; then
        echo "   ✓ $file → HTML generated"
    else
        echo "   ✗ $file → HTML NOT found"
    fi
done

# Check index includes new pages
echo -e "\n2. Checking index integration..."
if grep -q "assign_host.html" _build/html/index.html; then
    echo "   ✓ assign_host linked in index"
else
    echo "   ✗ assign_host NOT linked in index"
fi

if grep -q "quickstart_assign_host.html" _build/html/index.html; then
    echo "   ✓ quickstart_assign_host linked in index"
else
    echo "   ✗ quickstart_assign_host NOT linked in index"
fi

# Check cross-references
echo -e "\n3. Checking cross-references..."
if grep -q "assign_host.html" _build/html/simulations.html; then
    echo "   ✓ simulations.rst → assign_host.rst"
else
    echo "   ✗ simulations.rst ↛ assign_host.rst"
fi

if grep -q "simulations.html" _build/html/assign_host.html; then
    echo "   ✓ assign_host.rst → simulations.rst"
else
    echo "   ✗ assign_host.rst ↛ simulations.rst"
fi

# Check key sections
echo -e "\n4. Checking key sections in assign_host.html..."
sections=(
    "Overview"
    "Galaxy Catalog Requirements"
    "Algorithm Details"
    "API Reference"
    "assign_frbs_to_hosts"
)

for section in "${sections[@]}"; do
    if grep -q "$section" _build/html/assign_host.html; then
        echo "   ✓ $section"
    else
        echo "   ✗ $section missing"
    fi
done

# Check code examples
echo -e "\n5. Checking code examples..."
code_count=$(grep -c "highlight-python" _build/html/assign_host.html)
echo "   Found $code_count Python code blocks in assign_host.html"

if [ $code_count -gt 5 ]; then
    echo "   ✓ Sufficient code examples"
else
    echo "   ✗ Insufficient code examples"
fi

# Build status
echo -e "\n6. Build status..."
if [ -f "_build/html/assign_host.html" ] && [ -f "_build/html/quickstart_assign_host.html" ]; then
    echo "   ✓ Documentation built successfully"
else
    echo "   ✗ Documentation build incomplete"
fi

echo -e "\n=========================================="
echo "Verification complete!"
echo "=========================================="
