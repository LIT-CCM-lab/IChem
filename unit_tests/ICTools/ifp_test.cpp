#include <cassert>
#include <iostream>

extern unsigned chooseIFPType_forTests(bool polarOnly, bool extended, bool metalOnly);

// Check all 8 combinations
static void test_chooseIFPType_truth_table() {
    // 000: standard
    assert(chooseIFPType_forTests(false, false, false) == 0u);

    // 100: polar only => 1
    assert(chooseIFPType_forTests(true,  false, false) == 1u);

    // 010: extended only => 2
    assert(chooseIFPType_forTests(false, true,  false) == 2u);

    // 001: metal only => 4
    assert(chooseIFPType_forTests(false, false, true)  == 4u);

    // 110: polar + extended (no metal) => 3
    assert(chooseIFPType_forTests(true,  true,  false) == 3u);

    // 101: polar + metal, no extended => 1
    assert(chooseIFPType_forTests(true,  false, true)  == 1u);

    // 011: extended + metal, no polar => 2
    assert(chooseIFPType_forTests(false, true,  true)  == 2u);

    // 111: polar + extended + metal => 4
    assert(chooseIFPType_forTests(true,  true,  true)  == 4u);
}

int main() {
    test_chooseIFPType_truth_table();

    std::cout << "[OK] chooseIFPType tests passed\n";
    return 0;
}
