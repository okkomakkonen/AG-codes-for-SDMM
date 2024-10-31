from generate_plot import GASP, PoleGap

counts = {"total": 0, "better": 0, "equal": 0, "worse": 0}

K_limit = 50
X_limit = 50

# We test L <= K

for K in range(2, K_limit + 1, 2):
    for L in range(1, K + 1):
        for X in range(1, X_limit + 1):
            counts["total"] += 1

            R_PoleGap = PoleGap(K, L, X)
            R_GASP = GASP(K, L, X)

            if R_PoleGap < R_GASP:
                counts["better"] += 1
            elif R_PoleGap == R_GASP:
                counts["equal"] += 1
            elif R_PoleGap > R_GASP:
                counts["worse"] += 1
            else:
                raise

assert counts["better"] + counts["worse"] + counts["equal"] == counts["total"]

print(counts)
print(f"PoleGap is better than GASP in {100 * counts['better'] / counts['total']}\% of cases")
