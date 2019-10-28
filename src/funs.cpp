#include <exporter.h>
#include <RcppEigen.h>
// Code generated by Stan version 2.19.1

#include <stan/model/standalone_functions_header.hpp>

namespace user_f0820e0d7a9b64beafd46985e1dbe478 { 
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using namespace stan::math;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vector_d;
typedef Eigen::Matrix<double, 1, Eigen::Dynamic> row_vector_d;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matrix_d;

stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "unknown file name");
    reader.add_event(278, 276, "end", "unknown file name");
    return reader;
}

template <typename T0__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1>
lambdas(const T0__& L,
            const int& M, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 4;
        validate_non_negative_index("lambdas", "M", M);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> lambdas(M);
        stan::math::initialize(lambdas, DUMMY_VAR__);
        stan::math::fill(lambdas, DUMMY_VAR__);


        current_statement_begin__ = 5;
        for (int m = 1; m <= M; ++m) {

            current_statement_begin__ = 6;
            stan::model::assign(lambdas, 
                        stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                        pow(((stan::math::pi() * m) / (2 * L)), 2), 
                        "assigning variable lambdas");
        }
        current_statement_begin__ = 8;
        return stan::math::promote_scalar<fun_return_scalar_t__>(lambdas);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct lambdas_functor__ {
    template <typename T0__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1>
    operator()(const T0__& L,
            const int& M, std::ostream* pstream__) const {
        return lambdas(L, M, pstream__);
    }
};

template <typename T0__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T2__>::type, Eigen::Dynamic, 1>
basis_phi(const T0__& L,
              const int& m,
              const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& x, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 12;
        int N(0);
        (void) N;  // dummy to suppress unused var warning
        stan::math::fill(N, std::numeric_limits<int>::min());
        stan::math::assign(N,rows(x));

        current_statement_begin__ = 13;
        validate_non_negative_index("basis_phi", "N", N);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> basis_phi(N);
        stan::math::initialize(basis_phi, DUMMY_VAR__);
        stan::math::fill(basis_phi, DUMMY_VAR__);
        stan::math::assign(basis_phi,multiply((1 / stan::math::sqrt(L)), stan::math::sin(multiply(((stan::math::pi() * m) / (2 * L)), add(x, L)))));


        current_statement_begin__ = 14;
        return stan::math::promote_scalar<fun_return_scalar_t__>(basis_phi);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct basis_phi_functor__ {
    template <typename T0__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T2__>::type, Eigen::Dynamic, 1>
    operator()(const T0__& L,
              const int& m,
              const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& x, std::ostream* pstream__) const {
        return basis_phi(L, m, x, pstream__);
    }
};

template <typename T0__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T2__>::type, Eigen::Dynamic, Eigen::Dynamic>
basis_phis(const T0__& L,
               const int& M,
               const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& x, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 20;
        int N(0);
        (void) N;  // dummy to suppress unused var warning
        stan::math::fill(N, std::numeric_limits<int>::min());
        stan::math::assign(N,rows(x));

        current_statement_begin__ = 21;
        validate_non_negative_index("phis", "N", N);
        validate_non_negative_index("phis", "M", M);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> phis(N, M);
        stan::math::initialize(phis, DUMMY_VAR__);
        stan::math::fill(phis, DUMMY_VAR__);


        current_statement_begin__ = 22;
        for (int m = 1; m <= M; ++m) {

            current_statement_begin__ = 23;
            stan::model::assign(phis, 
                        stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list())), 
                        basis_phi(L, m, x, pstream__), 
                        "assigning variable phis");
        }
        current_statement_begin__ = 25;
        return stan::math::promote_scalar<fun_return_scalar_t__>(phis);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct basis_phis_functor__ {
    template <typename T0__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T2__>::type, Eigen::Dynamic, Eigen::Dynamic>
    operator()(const T0__& L,
               const int& M,
               const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& x, std::ostream* pstream__) const {
        return basis_phis(L, M, x, pstream__);
    }
};

template <typename T0__, typename T1__, typename T2__>
typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
spd(const T0__& alpha,
        const T1__& rho,
        const T2__& lambda, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 30;
        local_scalar_t__ spd(DUMMY_VAR__);
        (void) spd;  // dummy to suppress unused var warning
        stan::math::initialize(spd, DUMMY_VAR__);
        stan::math::fill(spd, DUMMY_VAR__);
        stan::math::assign(spd,(((pow(alpha, 2) * stan::math::sqrt((2 * stan::math::pi()))) * rho) * stan::math::exp(((-(.5) * pow(rho, 2)) * pow(lambda, 2)))));


        current_statement_begin__ = 31;
        return stan::math::promote_scalar<fun_return_scalar_t__>(spd);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct spd_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__>::type
    operator()(const T0__& alpha,
        const T1__& rho,
        const T2__& lambda, std::ostream* pstream__) const {
        return spd(alpha, rho, lambda, pstream__);
    }
};

template <typename T0__, typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
spds(const T0__& alpha,
         const T1__& rho,
         const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& lambdas, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 36;
        int M(0);
        (void) M;  // dummy to suppress unused var warning
        stan::math::fill(M, std::numeric_limits<int>::min());
        stan::math::assign(M,rows(lambdas));

        current_statement_begin__ = 37;
        validate_non_negative_index("spds", "M", M);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> spds(M);
        stan::math::initialize(spds, DUMMY_VAR__);
        stan::math::fill(spds, DUMMY_VAR__);


        current_statement_begin__ = 38;
        for (int m = 1; m <= M; ++m) {

            current_statement_begin__ = 39;
            stan::model::assign(spds, 
                        stan::model::cons_list(stan::model::index_uni(m), stan::model::nil_index_list()), 
                        spd(alpha, rho, get_base1(lambdas, m, "lambdas", 1), pstream__), 
                        "assigning variable spds");
        }
        current_statement_begin__ = 41;
        return stan::math::promote_scalar<fun_return_scalar_t__>(spds);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct spds_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
    operator()(const T0__& alpha,
         const T1__& rho,
         const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& lambdas, std::ostream* pstream__) const {
        return spds(alpha, rho, lambdas, pstream__);
    }
};

template <typename T0__, typename T1__, typename T2__, typename T3__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type, Eigen::Dynamic, Eigen::Dynamic>
spd_gp(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& x_phi,
           const T1__& alpha,
           const T2__& rho,
           const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& lambdas, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 47;
        int M(0);
        (void) M;  // dummy to suppress unused var warning
        stan::math::fill(M, std::numeric_limits<int>::min());
        stan::math::assign(M,cols(x_phi));

        current_statement_begin__ = 48;
        validate_non_negative_index("spdsOut", "M", M);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> spdsOut(M);
        stan::math::initialize(spdsOut, DUMMY_VAR__);
        stan::math::fill(spdsOut, DUMMY_VAR__);
        stan::math::assign(spdsOut,spds(alpha, rho, lambdas, pstream__));


        current_statement_begin__ = 49;
        return stan::math::promote_scalar<fun_return_scalar_t__>(diag_post_multiply(x_phi, stan::math::sqrt(spdsOut)));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct spd_gp_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type, Eigen::Dynamic, Eigen::Dynamic>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& x_phi,
           const T1__& alpha,
           const T2__& rho,
           const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& lambdas, std::ostream* pstream__) const {
        return spd_gp(x_phi, alpha, rho, lambdas, pstream__);
    }
};

template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type, Eigen::Dynamic, 1>
spd_gp_fast(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& x_phi,
                const T1__& alpha,
                const T2__& rho,
                const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& lambdas,
                const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& gp_z, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 54;
        int M(0);
        (void) M;  // dummy to suppress unused var warning
        stan::math::fill(M, std::numeric_limits<int>::min());
        stan::math::assign(M,cols(x_phi));

        current_statement_begin__ = 55;
        validate_non_negative_index("spdsOut", "M", M);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> spdsOut(M);
        stan::math::initialize(spdsOut, DUMMY_VAR__);
        stan::math::fill(spdsOut, DUMMY_VAR__);
        stan::math::assign(spdsOut,spds(alpha, rho, lambdas, pstream__));


        current_statement_begin__ = 56;
        return stan::math::promote_scalar<fun_return_scalar_t__>(multiply(x_phi, elt_multiply(spdsOut, gp_z)));
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct spd_gp_fast_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__, typename T4__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__, typename boost::math::tools::promote_args<T4__>::type>::type, Eigen::Dynamic, 1>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& x_phi,
                const T1__& alpha,
                const T2__& rho,
                const Eigen::Matrix<T3__, Eigen::Dynamic, 1>& lambdas,
                const Eigen::Matrix<T4__, Eigen::Dynamic, 1>& gp_z, std::ostream* pstream__) const {
        return spd_gp_fast(x_phi, alpha, rho, lambdas, gp_z, pstream__);
    }
};

template <typename T0__, typename T3__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T3__>::type, Eigen::Dynamic, Eigen::Dynamic>
omega_one(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& lambda_loc_mat,
              const std::vector<std::vector<int> >& F_inds,
              const std::vector<int>& F_inds_num,
              const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& shat, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 61;
        int N(0);
        (void) N;  // dummy to suppress unused var warning
        stan::math::fill(N, std::numeric_limits<int>::min());
        stan::math::assign(N,rows(shat));

        current_statement_begin__ = 62;
        int F(0);
        (void) F;  // dummy to suppress unused var warning
        stan::math::fill(F, std::numeric_limits<int>::min());
        stan::math::assign(F,rows(lambda_loc_mat));

        current_statement_begin__ = 63;
        int J(0);
        (void) J;  // dummy to suppress unused var warning
        stan::math::fill(J, std::numeric_limits<int>::min());
        stan::math::assign(J,cols(shat));

        current_statement_begin__ = 64;
        validate_non_negative_index("vhat", "N", N);
        validate_non_negative_index("vhat", "J", J);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> vhat(N, J);
        stan::math::initialize(vhat, DUMMY_VAR__);
        stan::math::fill(vhat, DUMMY_VAR__);
        stan::math::assign(vhat,elt_multiply(shat, shat));

        current_statement_begin__ = 65;
        validate_non_negative_index("vhat_sum", "N", N);
        validate_non_negative_index("vhat_sum", "F", F);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> vhat_sum(N, F);
        stan::math::initialize(vhat_sum, DUMMY_VAR__);
        stan::math::fill(vhat_sum, DUMMY_VAR__);

        current_statement_begin__ = 66;
        validate_non_negative_index("numerator", "F", F);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> numerator(F);
        stan::math::initialize(numerator, DUMMY_VAR__);
        stan::math::fill(numerator, DUMMY_VAR__);

        current_statement_begin__ = 67;
        validate_non_negative_index("omega", "N", N);
        validate_non_negative_index("omega", "F", F);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> omega(N, F);
        stan::math::initialize(omega, DUMMY_VAR__);
        stan::math::fill(omega, DUMMY_VAR__);


        current_statement_begin__ = 68;
        for (int f = 1; f <= F; ++f) {

            current_statement_begin__ = 69;
            stan::model::assign(vhat_sum, 
                        stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(f), stan::model::nil_index_list())), 
                        multiply(stan::model::rvalue(vhat, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(F_inds, stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_min_max(1, get_base1(F_inds_num, f, "F_inds_num", 1)), stan::model::nil_index_list())), "F_inds")), stan::model::nil_index_list())), "vhat"), rep_vector(1, get_base1(F_inds_num, f, "F_inds_num", 1))), 
                        "assigning variable vhat_sum");
            current_statement_begin__ = 70;
            stan::model::assign(numerator, 
                        stan::model::cons_list(stan::model::index_uni(f), stan::model::nil_index_list()), 
                        pow(sum(stan::model::rvalue(lambda_loc_mat, stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "lambda_loc_mat")), 2), 
                        "assigning variable numerator");
            current_statement_begin__ = 71;
            stan::model::assign(omega, 
                        stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(f), stan::model::nil_index_list())), 
                        elt_divide(get_base1(numerator, f, "numerator", 1), add(get_base1(numerator, f, "numerator", 1), stan::model::rvalue(vhat_sum, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(f), stan::model::nil_index_list())), "vhat_sum"))), 
                        "assigning variable omega");
        }
        current_statement_begin__ = 73;
        return stan::math::promote_scalar<fun_return_scalar_t__>(omega);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct omega_one_functor__ {
    template <typename T0__, typename T3__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T3__>::type, Eigen::Dynamic, Eigen::Dynamic>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& lambda_loc_mat,
              const std::vector<std::vector<int> >& F_inds,
              const std::vector<int>& F_inds_num,
              const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& shat, std::ostream* pstream__) const {
        return omega_one(lambda_loc_mat, F_inds, F_inds_num, shat, pstream__);
    }
};

Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
loadings_to_ones(const std::vector<std::vector<int> >& F_inds,
                     const std::vector<int>& F_inds_num, std::ostream* pstream__) {
    typedef double local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 79;
        int F(0);
        (void) F;  // dummy to suppress unused var warning
        stan::math::fill(F, std::numeric_limits<int>::min());
        stan::math::assign(F,size(F_inds));

        current_statement_begin__ = 80;
        int J(0);
        (void) J;  // dummy to suppress unused var warning
        stan::math::fill(J, std::numeric_limits<int>::min());
        stan::math::assign(J,size(get_base1(F_inds, 1, "F_inds", 1)));

        current_statement_begin__ = 81;
        validate_non_negative_index("lambda_ones", "F", F);
        validate_non_negative_index("lambda_ones", "J", J);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> lambda_ones(F, J);
        stan::math::initialize(lambda_ones, DUMMY_VAR__);
        stan::math::fill(lambda_ones, DUMMY_VAR__);


        current_statement_begin__ = 82;
        for (int f = 1; f <= F; ++f) {

            current_statement_begin__ = 83;
            for (int j = 1; j <= J; ++j) {

                current_statement_begin__ = 84;
                stan::model::assign(lambda_ones, 
                            stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_uni(j), stan::model::nil_index_list())), 
                            0, 
                            "assigning variable lambda_ones");
            }
            current_statement_begin__ = 86;
            stan::model::assign(lambda_ones, 
                        stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(F_inds, stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_min_max(1, get_base1(F_inds_num, f, "F_inds_num", 1)), stan::model::nil_index_list())), "F_inds")), stan::model::nil_index_list())), 
                        rep_row_vector(1, get_base1(F_inds_num, f, "F_inds_num", 1)), 
                        "assigning variable lambda_ones");
        }
        current_statement_begin__ = 88;
        return stan::math::promote_scalar<fun_return_scalar_t__>(lambda_ones);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct loadings_to_ones_functor__ {
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    operator()(const std::vector<std::vector<int> >& F_inds,
                     const std::vector<int>& F_inds_num, std::ostream* pstream__) const {
        return loadings_to_ones(F_inds, F_inds_num, pstream__);
    }
};

Eigen::Matrix<double, Eigen::Dynamic, 1>
ones(const int& num, std::ostream* pstream__) {
    typedef double local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 92;
        validate_non_negative_index("ones", "num", num);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> ones(num);
        stan::math::initialize(ones, DUMMY_VAR__);
        stan::math::fill(ones, DUMMY_VAR__);


        current_statement_begin__ = 93;
        for (int n = 1; n <= num; ++n) {

            current_statement_begin__ = 94;
            stan::model::assign(ones, 
                        stan::model::cons_list(stan::model::index_uni(n), stan::model::nil_index_list()), 
                        1, 
                        "assigning variable ones");
        }
        current_statement_begin__ = 96;
        return stan::math::promote_scalar<fun_return_scalar_t__>(ones);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct ones_functor__ {
            Eigen::Matrix<double, Eigen::Dynamic, 1>
    operator()(const int& num, std::ostream* pstream__) const {
        return ones(num, pstream__);
    }
};

template <typename T0__, typename T3__, typename T4__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T3__, T4__>::type, Eigen::Dynamic, Eigen::Dynamic>
omega_two(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& lambda_loc_mat,
              const std::vector<std::vector<int> >& F_inds,
              const std::vector<int>& F_inds_num,
              const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& theta_cor_L,
              const Eigen::Matrix<T4__, Eigen::Dynamic, Eigen::Dynamic>& shat, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T3__, T4__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning

    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 101;
        int N(0);
        (void) N;  // dummy to suppress unused var warning
        stan::math::fill(N, std::numeric_limits<int>::min());
        stan::math::assign(N,rows(shat));

        current_statement_begin__ = 102;
        int F(0);
        (void) F;  // dummy to suppress unused var warning
        stan::math::fill(F, std::numeric_limits<int>::min());
        stan::math::assign(F,rows(lambda_loc_mat));

        current_statement_begin__ = 103;
        int J(0);
        (void) J;  // dummy to suppress unused var warning
        stan::math::fill(J, std::numeric_limits<int>::min());
        stan::math::assign(J,cols(shat));

        current_statement_begin__ = 104;
        validate_non_negative_index("vhat", "N", N);
        validate_non_negative_index("vhat", "J", J);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> vhat(N, J);
        stan::math::initialize(vhat, DUMMY_VAR__);
        stan::math::fill(vhat, DUMMY_VAR__);
        stan::math::assign(vhat,elt_multiply(shat, shat));

        current_statement_begin__ = 105;
        validate_non_negative_index("theta_loc_cov", "F", F);
        validate_non_negative_index("theta_loc_cov", "F", F);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> theta_loc_cov(F, F);
        stan::math::initialize(theta_loc_cov, DUMMY_VAR__);
        stan::math::fill(theta_loc_cov, DUMMY_VAR__);
        stan::math::assign(theta_loc_cov,multiply_lower_tri_self_transpose(theta_cor_L));

        current_statement_begin__ = 108;
        validate_non_negative_index("implied_cov_fixed", "J", J);
        validate_non_negative_index("implied_cov_fixed", "J", J);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> implied_cov_fixed(J, J);
        stan::math::initialize(implied_cov_fixed, DUMMY_VAR__);
        stan::math::fill(implied_cov_fixed, DUMMY_VAR__);
        stan::math::assign(implied_cov_fixed,multiply(multiply(transpose(lambda_loc_mat), theta_loc_cov), lambda_loc_mat));

        current_statement_begin__ = 109;
        validate_non_negative_index("lambda_ones", "F", F);
        validate_non_negative_index("lambda_ones", "J", J);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> lambda_ones(F, J);
        stan::math::initialize(lambda_ones, DUMMY_VAR__);
        stan::math::fill(lambda_ones, DUMMY_VAR__);
        stan::math::assign(lambda_ones,loadings_to_ones(F_inds, F_inds_num, pstream__));

        current_statement_begin__ = 110;
        validate_non_negative_index("numerator", "F", F);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> numerator(F);
        stan::math::initialize(numerator, DUMMY_VAR__);
        stan::math::fill(numerator, DUMMY_VAR__);

        current_statement_begin__ = 111;
        validate_non_negative_index("omega", "N", N);
        validate_non_negative_index("omega", "F", F);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, Eigen::Dynamic> omega(N, F);
        stan::math::initialize(omega, DUMMY_VAR__);
        stan::math::fill(omega, DUMMY_VAR__);


        current_statement_begin__ = 112;
        for (int f = 1; f <= F; ++f) {

            current_statement_begin__ = 113;
            stan::model::assign(numerator, 
                        stan::model::cons_list(stan::model::index_uni(f), stan::model::nil_index_list()), 
                        pow(sum(stan::model::rvalue(lambda_loc_mat, stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_omni(), stan::model::nil_index_list())), "lambda_loc_mat")), 2), 
                        "assigning variable numerator");
            current_statement_begin__ = 114;
            stan::model::assign(omega, 
                        stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_uni(f), stan::model::nil_index_list())), 
                        elt_divide(get_base1(numerator, f, "numerator", 1), add(multiply(multiply(get_base1(lambda_ones, f, "lambda_ones", 1), implied_cov_fixed), transpose(get_base1(lambda_ones, f, "lambda_ones", 1))), multiply(stan::model::rvalue(vhat, stan::model::cons_list(stan::model::index_omni(), stan::model::cons_list(stan::model::index_multi(stan::model::rvalue(F_inds, stan::model::cons_list(stan::model::index_uni(f), stan::model::cons_list(stan::model::index_min_max(1, get_base1(F_inds_num, f, "F_inds_num", 1)), stan::model::nil_index_list())), "F_inds")), stan::model::nil_index_list())), "vhat"), ones(get_base1(F_inds_num, f, "F_inds_num", 1), pstream__)))), 
                        "assigning variable omega");
        }
        current_statement_begin__ = 116;
        return stan::math::promote_scalar<fun_return_scalar_t__>(omega);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}


struct omega_two_functor__ {
    template <typename T0__, typename T3__, typename T4__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T3__, T4__>::type, Eigen::Dynamic, Eigen::Dynamic>
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, Eigen::Dynamic>& lambda_loc_mat,
              const std::vector<std::vector<int> >& F_inds,
              const std::vector<int>& F_inds_num,
              const Eigen::Matrix<T3__, Eigen::Dynamic, Eigen::Dynamic>& theta_cor_L,
              const Eigen::Matrix<T4__, Eigen::Dynamic, Eigen::Dynamic>& shat, std::ostream* pstream__) const {
        return omega_two(lambda_loc_mat, F_inds, F_inds_num, theta_cor_L, shat, pstream__);
    }
};

 } 
// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1>
lambdas(const double& L,
            const int& M){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::lambdas<double>(L, M, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1>
basis_phi(const double& L,
              const int& m,
              const Eigen::Matrix<double, Eigen::Dynamic, 1>& x){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::basis_phi<double, double>(L, m, x, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
basis_phis(const double& L,
               const int& M,
               const Eigen::Matrix<double, Eigen::Dynamic, 1>& x){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::basis_phis<double, double>(L, M, x, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
double
spd(const double& alpha,
        const double& rho,
        const double& lambda){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::spd<double, double, double>(alpha, rho, lambda, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1>
spds(const double& alpha,
         const double& rho,
         const Eigen::Matrix<double, Eigen::Dynamic, 1>& lambdas){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::spds<double, double, double>(alpha, rho, lambdas, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
spd_gp(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& x_phi,
           const double& alpha,
           const double& rho,
           const Eigen::Matrix<double, Eigen::Dynamic, 1>& lambdas){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::spd_gp<double, double, double, double>(x_phi, alpha, rho, lambdas, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1>
spd_gp_fast(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& x_phi,
                const double& alpha,
                const double& rho,
                const Eigen::Matrix<double, Eigen::Dynamic, 1>& lambdas,
                const Eigen::Matrix<double, Eigen::Dynamic, 1>& gp_z){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::spd_gp_fast<double, double, double, double, double>(x_phi, alpha, rho, lambdas, gp_z, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
omega_one(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& lambda_loc_mat,
              const std::vector<std::vector<int> >& F_inds,
              const std::vector<int>& F_inds_num,
              const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& shat){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::omega_one<double, double>(lambda_loc_mat, F_inds, F_inds_num, shat, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
loadings_to_ones(const std::vector<std::vector<int> >& F_inds,
                     const std::vector<int>& F_inds_num){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::loadings_to_ones(F_inds, F_inds_num, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, 1>
ones(const int& num){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::ones(num, &Rcpp::Rcout);
}

// [[Rcpp::depends(rstan)]]
// [[Rcpp::export]]
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
omega_two(const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& lambda_loc_mat,
              const std::vector<std::vector<int> >& F_inds,
              const std::vector<int>& F_inds_num,
              const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& theta_cor_L,
              const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& shat){
  return 
user_f0820e0d7a9b64beafd46985e1dbe478::omega_two<double, double, double>(lambda_loc_mat, F_inds, F_inds_num, theta_cor_L, shat, &Rcpp::Rcout);
}

