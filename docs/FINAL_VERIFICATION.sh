#!/bin/bash
echo "================================================================================"
echo "                    FINAL VERIFICATION CHECKLIST"
echo "================================================================================"
echo

# Check module
echo "✓ Module Implementation"
if [ -f "astropath/simulations/assign_host.py" ]; then
    lines=$(wc -l < astropath/simulations/assign_host.py)
    echo "  - assign_host.py: $lines lines"
else
    echo "  ✗ assign_host.py NOT FOUND"
fi

# Check __init__.py
echo
echo "✓ Module Integration"
if grep -q "assign_frbs_to_hosts" astropath/simulations/__init__.py; then
    echo "  - Functions exported in __init__.py"
else
    echo "  ✗ Functions NOT exported"
fi

# Check RST docs
echo
echo "✓ RST Documentation"
for file in assign_host.rst quickstart_assign_host.rst; do
    if [ -f "docs/$file" ]; then
        size=$(ls -lh "docs/$file" | awk '{print $5}')
        echo "  - $file: $size"
    else
        echo "  ✗ $file NOT FOUND"
    fi
done

# Check notebook
echo
echo "✓ Jupyter Notebook"
if [ -f "docs/nb/Simulate_Assign_FRBs.ipynb" ]; then
    cells=$(grep -c '"cell_type"' docs/nb/Simulate_Assign_FRBs.ipynb)
    echo "  - Simulate_Assign_FRBs.ipynb: $cells cells"
else
    echo "  ✗ Notebook NOT FOUND"
fi

# Check tests
echo
echo "✓ Test Files"
if [ -f "test_assign_host.py" ]; then
    echo "  - test_assign_host.py: present"
else
    echo "  ✗ test_assign_host.py NOT FOUND"
fi
if [ -f "examples/assign_frbs_example.py" ]; then
    echo "  - examples/assign_frbs_example.py: present"
else
    echo "  ✗ examples/assign_frbs_example.py NOT FOUND"
fi

# Check index integration
echo
echo "✓ Documentation Index"
if grep -q "assign_host" docs/index.rst; then
    echo "  - assign_host in index"
fi
if grep -q "quickstart_assign_host" docs/index.rst; then
    echo "  - quickstart_assign_host in index"
fi
if grep -q "Simulate_Assign_FRBs" docs/index.rst; then
    echo "  - Simulate_Assign_FRBs in index"
fi

# Check cross-references
echo
echo "✓ Cross-References"
if grep -q "assign_host" docs/simulations.rst; then
    echo "  - simulations.rst → assign_host.rst"
fi

# Check HTML build
echo
echo "✓ HTML Documentation"
if [ -f "_build/html/assign_host.html" ]; then
    size=$(ls -lh "_build/html/assign_host.html" | awk '{print $5}')
    echo "  - assign_host.html: $size"
fi
if [ -f "_build/html/quickstart_assign_host.html" ]; then
    size=$(ls -lh "_build/html/quickstart_assign_host.html" | awk '{print $5}')
    echo "  - quickstart_assign_host.html: $size"
fi

echo
echo "================================================================================"
echo "                         VERIFICATION COMPLETE"
echo "================================================================================"
