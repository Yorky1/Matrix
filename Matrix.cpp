#include <vector>
#include <iostream>
#include <exception>

class Matrix {
private:
    std::vector<std::vector<Rational>> data;
    int N, M;

public:
    Matrix(int n, int m, const std::vector<std::vector<Rational>>& a) : N(n), M(m), data(a) {}

    Matrix(int n = 0, int m = 0) : N(n), M(m) {
        data.resize(n, std::vector<Rational>(m));
    }

    void resize(int n, int m) {
        N = n;
        M = m;
        data.resize(n, std::vector<Rational>(m));
    }

    std::pair<int, int> size() const {
        return { N, M };
    }

    Matrix operator*(const Matrix& other) {
        if (M != other.N) {
            return {};
        }
        Matrix result(N, other.M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < other.M; ++j) {
                for (int k = 0; k < M; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }

    Rational trace() {
        Rational result = 0;
        for (int i = 0; i < min(N, M); ++i) {
            result += data[i][i];
        }
        return result;
    }

    Matrix operator+ (const Matrix& other) const {
        Matrix result(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                result.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return result;
    }

    Matrix operator- (const Matrix& other) const {
        Matrix result(N, M);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                result.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return result;
    }

    Matrix T() {
        Matrix result(M, N);
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    std::vector<Rational>& operator[] (int i) {
        return data[i];
    }
};


std::ostream& operator<<(std::ostream& out, Matrix& mat) {
    auto [n, m] = mat.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            out << mat[i][j] << " ";
        }
        out << '\n';
    }
    return out;
}

std::istream& operator>>(std::istream& in, Matrix& mat) {
    int n, m;
    in >> n >> m;
    mat.resize(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cin >> mat[i][j];
        }
    }
    return in;
}

Matrix ones(int n) {
    Matrix res(n, n);
    for (int i = 0; i < n; ++i) {
        res[i][i] = 1;
    }
    return res;
}

Matrix getReversed(Matrix a) {
    auto [n, m] = a.size();

    if (n != m) {
        throw exception("There is no inverse matrix");
    }

    Matrix res = ones(n);

    for (int col = 0, row = 0; col < m && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i)
            if (abs(a[i][col]) > abs(a[sel][col]))
                sel = i;
        if (a[sel][col] == 0)
            continue;

        swap(a[sel], a[row]);
        swap(res[sel], res[row]);

        for (int i = 0; i < n; ++i) {
            if (i != row) {
                Rational c = a[i][col] / a[row][col];
                for (int j = col; j < m; ++j)
                    a[i][j] -= a[row][j] * c;
                for (int j = 0; j < m; ++j) {
                    res[i][j] -= res[row][j] * c;
                }
            }
        }
        ++row;
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            res[i][j] /= a[i][i];
        }
        a[i][i] /= a[i][i];
    }

    for (int i = 0; i < n; ++i) {
        if (a[i][i] != 1) {
            throw exception("There is no inverse matrix");
        }
    }

    return res;
}

Matrix getEchelon(Matrix a) {
    auto [n, m] = a.size();

    for (int col = 0, row = 0; col < m && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i)
            if (abs(a[i][col]) > abs(a[sel][col]))
                sel = i;
        if (a[sel][col] == 0)
            continue;

        swap(a[sel], a[row]);

        for (int i = row + 1; i < n; ++i) {
            Rational c = a[i][col] / a[row][col];
            for (int j = col; j < m; ++j)
                a[i][j] -= a[row][j] * c;
        }
        ++row;
    }

    return a;
}

Matrix getReduced(Matrix a) {
    auto [n, m] = a.size();

    for (int col = 0, row = 0; col < m && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i)
            if (abs(a[i][col]) > abs(a[sel][col]))
                sel = i;
        if (a[sel][col] == 0)
            continue;

        swap(a[sel], a[row]);

        {
            Rational c = a[row][col];
            for (int i = col; i < m; ++i) {
                a[row][i] /= c;
            }
        }

        for (int i = 0; i < n; ++i) {
            if (i == row) continue;
            Rational c = a[i][col] / a[row][col];
            for (int j = col; j < m; ++j)
                a[i][j] -= a[row][j] * c;
        }
        ++row;
    }

    return a;
}
