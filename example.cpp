#include "pohlig_hellman.h"
#include <iostream>
#include <chrono>

static std::chrono::steady_clock::time_point t_start;
static void timer_start() {
    t_start = std::chrono::steady_clock::now();
}
static void timer_stop() {
    auto t_end = std::chrono::steady_clock::now();
    int elapsed_s = (int)std::chrono::duration_cast<std::chrono::seconds>(t_end - t_start).count();
    int elapsed_ms = (int)std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();
    elapsed_ms -= elapsed_s * 1000;

    std::cout << "Took " << elapsed_s << "s " << elapsed_ms << "ms\n";
}

void test_xoro_distance(uint64_t target_dist) {
    timer_start();
    std::cout << "----- Xoro distance test\n";

    Xoroshiro a = {123456789, 0};
    Xoroshiro b = a;

    FXTMatrix adv, advT;
    fastXoroMatrixPower(&XOROSHIRO_STANDARD_FXTM, target_dist, &adv);
    transposeFXTM(&adv, &advT);
    fastAdvanceXoroshiroFXTM(&b, &advT);

    InfInt dist = get_xoroshiro_state_distance(a, b, true);
    std::cout << "Result:\ndistance = " << dist << "\n";
    timer_stop();
    std::cout << '\n';
}

void test_matrix_solver() {
    timer_start();
    std::cout << "----- Matrix column solver test\n";

    Xoroshiro target_col = {0, 1};
    InfInt dist = get_xoroshiro_jump_with_column(0, target_col, true);

    std::cout << "Result:\ndistance = " << dist << "\n";

    bool any_different = false;

    Xoroshiro a = {987654321, 123456789};
    for (int i = 0; i < 256; i++) {
        Xoroshiro b = a;
        advance_infint(&b, dist);

        if ((b.hi & 1) != (a.lo & 1)) {
            any_different = true;
        }

        xNextLong(&a);
    }

    std::cout << "LSB of state.lo = LSB of adv(state, distance).hi ?  " 
            << (any_different ? "false" : "true") << "\n";
    timer_stop();
    std::cout << '\n';
}

int main() {
    test_xoro_distance(444444444444444444);
    test_matrix_solver();
    return 0;
}