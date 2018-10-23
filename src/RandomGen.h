#ifndef WRS_RANDOMGEN_H
#define WRS_RANDOMGEN_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>

namespace rws {

// random device
    static std::random_device g_rd;

    enum class SamplingMethod {
        RouletteWheel,
        BiasedCoin,
        AliasMethod
    };

// default choosing Alias Method
    template<typename T, SamplingMethod algo = SamplingMethod::AliasMethod>
    class RandomGen;

    template<typename T, SamplingMethod algo>
    class RandomGenBase {
    public:
        // default choosing algo AliasMethod
        RandomGenBase(std::vector <T> const &vals, std::vector <T> const &probs)
                :
                m_mt{g_rd()}, m_dist(T(0), T(1)), m_vals{vals}, m_probs{probs} {
            init();
        }

    public:
        T next() {
            return gen();
        }

    private:
        virtual void init() {};

        virtual T gen() = 0;

    protected:
        void checkValidDistribution() {
            if (m_probs.size() == 0 || m_vals.size() == 0 || m_probs.size() != m_vals.size()) {
                throw std::runtime_error("Error : Invalid Values and Probabilities.");
            }

            if (std::abs(std::accumulate(m_probs.cbegin(), m_probs.cend(), T(0)) - T(1))
                > T(100) * std::numeric_limits<T>::epsilon()) {
                // choosing 'numeric_limits' for the error allowance
                throw std::runtime_error("Error : Probabilities should sum up to 1.");
            }
        }

    protected:
        // common data members shared by all algorithms
        std::vector <T> m_vals;
        std::vector <T> m_probs;

    protected:
        // use internal uniform distribution
        std::mt19937 m_mt;
        std::uniform_real_distribution <T> m_dist;
    };

/*
 * RouletteWheel
 */
    template<typename T>
    class RandomGen<T, SamplingMethod::RouletteWheel> : public RandomGenBase<T, SamplingMethod::RouletteWheel> {
    public:
        RandomGen(std::vector <T> const &vals, std::vector <T> const &probs)
                : RandomGenBase<T, SamplingMethod::RouletteWheel>(vals, probs), n(this->m_probs.size()), A(n) {
            init();
        }

    private:
        void init() override {
            this->checkValidDistribution();

            // calculate the cummulative distribution
            A[0] = this->m_probs[0];
            for (auto i = 1; i <= n - 1; i++) {
                A[i] = A[i - 1] + this->m_probs[i];
            }
        }

        T gen() override {
            auto x = this->m_dist(this->m_mt);
            // std::upper_bound is implemented by binary search
            auto res = std::upper_bound(A.cbegin(), A.cend(), x);
            return res == A.cbegin() ? this->m_vals[0] : this->m_vals.at(
                    std::distance(A.cbegin(), res));
        }

    private:
        std::size_t n;
        std::vector <T> A;
    };

/*
 * BiasedCoin
 */
    template<typename T>
    class RandomGen<T, SamplingMethod::BiasedCoin> : public RandomGenBase<T, SamplingMethod::BiasedCoin> {
    public:
        RandomGen(std::vector <T> const &vals, std::vector <T> const &probs)
                : RandomGenBase<T, SamplingMethod::BiasedCoin>(vals, probs) {
            init();
        }

    private:
        void init() override {
            this->checkValidDistribution();
        }

        T gen() override {
            auto mass = T{1};
            auto n = this->m_probs.size();

            for (auto i = 0; i <= n - 1; i++) {
                auto x = this->m_dist(this->m_mt);
                if (x < (this->m_probs[i] / mass)) {
                    return this->m_vals.at(i);
                } else {
                    mass -= this->m_probs[i];
                }
            }
            throw std::runtime_error("Error in Biased Coin sampling.");
        }
    };

/*
 * AliasMethod
 */
    template<typename T>
    class RandomGen<T, SamplingMethod::AliasMethod> : public RandomGenBase<T, SamplingMethod::AliasMethod> {
    public:
        RandomGen(std::vector <T> const &vals, std::vector <T> const &probs)
                : RandomGenBase<T, SamplingMethod::AliasMethod>(vals, probs), n{this->m_probs.size()}, alias(n),
                  prob(n) {
            init();
        }

    private:
        void init() override {
            this->checkValidDistribution();

            // helper stacks for calculating the alias
            std::stack <T> small;
            std::stack <T> large;

            // make sure 0 and 1 are correctly handled by Type.
            auto zero = T(0);
            auto one = T(1);

            // normalize probabilities
            auto scale = n;
            std::for_each(this->m_probs.begin(), this->m_probs.end(), [scale](T &e) {
                e *= scale;
            });

            // initialize small and large
            for (auto i = 0; i < n; i++) {
                if (this->m_probs[i] < one) {
                    small.push(i);
                } else {
                    large.push(i);
                }
            }

            // initialize alias list
            while (!small.empty() && !large.empty()) {
                auto l = small.top();
                small.pop();

                auto g = large.top();
                large.pop();

                prob[l] = this->m_probs[l];
                alias[l] = g;

                this->m_probs[g] = (this->m_probs[g] + this->m_probs[l]) - one;

                if (this->m_probs[g] < one) {
                    small.push(g);
                } else {
                    large.push(g);
                }
            }

            // edge case - due to numerical instability
            while (!large.empty()) {
                auto g = large.top();
                large.pop();

                prob[g] = one;
            }

            while (!small.empty()) {
                auto l = small.top();
                small.pop();

                prob[l] = one;
            }
        }

        T gen() override {
            auto x = this->m_dist(this->m_mt);
            auto i = std::size_t(std::floor(x * n));

            if (this->m_dist(this->m_mt) < prob[i]) {
                return this->m_vals[i];
            } else {
                return this->m_vals[alias[i]];
            }
        }

    private:
        std::size_t n;
        std::vector <std::size_t> alias;
        std::vector <T> prob;
    };

}

#endif //WRS_RANDOMGEN_H
