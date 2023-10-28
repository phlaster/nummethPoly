int N = 4; double delta = 1e-2;
non_singular_needed:
    Mtr B = randMtr(N);
if (det(B) < 1e-7)
    goto non_singular_needed; // det(B) > 0

Vec eigens = randVec(N);
Mtr D = diag(eigens);
