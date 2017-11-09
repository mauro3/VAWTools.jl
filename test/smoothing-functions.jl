@test max_smooth(-0.1, 100, 0.0001) ≈ 100
@test max_smooth(-0.1, -100, 0.0001) ≈ -0.1
@test max_smooth(-0.1, 100, 0) == max(-0.1,100)

@test min_smooth(-0.1, 100, 0.0001) ≈ -0.1
@test min_smooth(-0.1, -100, 0.0001) ≈ -100
