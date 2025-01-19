#include <boost/math/special_functions/erf.hpp>
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <iomanip>

#define M_PI 3.14159265358979323846
using namespace std;
using namespace boost::math;

int option;
char text[8];

// long double Combinationfact(int n, int k) {
//     if (k == 0 || k == n) {
//         return 1;
//     }
//     return Combinationfact(n - 1, k - 1) + Combinationfact(n - 1, k);
//     return n * Combinationfact(n - 1, 1) / (Combinationfact(stop, 1) * Combinationfact(n, n - stop));
// }

long double Combinationfact(int n, int k) {
    long double numerator = 1, denominator = 1;

    for (int i = 1; i <= k; ++i) {
        numerator *= n--;
        denominator *= i;
    }
    return numerator / denominator;
}

long double factorial(int n) {
    if (n == 0) {
        return 1;
    } return n * factorial(n - 1);
}

void BernouliFunction(float probability, int i) {
    double success, mean = 0.0, variance = 0.0;

    do {
        cout << "1. Probability Mass Function" << endl;
        cout << "2. Mean and Variance Function" << endl;
        cout << "Choose what to find: ";
        cin >> option;
        switch (option) {
        case 1:
            // Probability Mass Function
            if (i != 0) {
                success = pow(probability, 1 - i)*pow((1 - probability), i);
                cout << "Probability: " << success << endl;
            }
            break;
        case 2:
            // Mean and Variances Bernouli's Distribution
            if (!(mean) && !(variance)) {
                for (int i = 0; i <= 1; i++) {
                    mean += i * (probability);
                    variance += (pow(i, 2) * probability) - pow(mean, 2);
                }
            }

            cout << "Mean: " << mean << endl;
            cout << "Variance: " << variance << endl;
            mean = 0.0, variance = 0.0;
            break;
        default:
            break;
        }
        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
   } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
}

void BinomialFunction(float probability) {
    int n, k;
    long double BinomialForm = 0.0, mean = 0.0, variance = 0.0, CDF_Result = 0.0;

    do {
        cout << "1. Probability Mass Function" << endl;
        cout << "2. Mean and Variance Function" << endl;
        cout << "3. Cumulative Distributive Function" << endl;
        cout << "4. Two Cumulative Distributive Function" << endl;
        cout << "Choose what to find: ";
        cin >> option;
        switch (option){
        case 1:
            // Probaility Mass Function
            cout << "Total Object/Attempt/Experiment: ";
            cin >> n;
            cout << "Request Object To Find Probability: ";
            cin >> k;

            BinomialForm = Combinationfact(n, k) * pow(probability, k) * pow(1 - probability, n - k);
            cout << "PMF Result: " << BinomialForm << endl;
            break;

        case 2:
            // Mean and Variances Binomial Distribution
            cout << "Total Object/Attempt/Experiment: ";
            cin >> n;

            if (!(mean) && !(variance)) {
                mean = n * (probability);
                variance = mean * (1 - probability);
            }
            
            cout << "Mean: " << mean << endl;
            cout << "Variance: " << variance << endl;
            mean = 0.0, variance = 0.0;
            break;

        case 3:
            // Cumulative Distributive Function Binomial Distribution
            cout << "Total Object/Attempt/Experiment: ";
            cin >> n;
            cout << "Request Object To Find Probability: ";
            cin >> k;

            CDF_Result = 0.0;
            for (int i = 0; i <= k; i++){
                CDF_Result += Combinationfact(n, i) * pow(probability, i) * pow(1 - probability, n - i);
            }

            cout << "Result (X <= x): " << CDF_Result << endl;
            cout << "Result (X >= x): " << 1 - CDF_Result << endl;
            break;

        case 4:
            // Two Cumulative Distributive Function Binomial Distribution
            int k[2];
            cout << "Total Object/Attempt/Experiment: ";
            cin >> n;
            
            for (unsigned int i = 0; i < sizeof(k)/sizeof(int);i++) {
                cout << "P(X) k" << i + 1 << ": ";
                cin >> k[i];
            }

            if (k[0] > k[1]) {
                swap(k[0], k[1]);
            }
            
            for (int i = 0; i <= k[1]; i++) {
                double currentP = Combinationfact(n, i) * pow(probability, i) * pow(1 - probability, n - i);
                if (i <= 1 - k[0]) {
                    BinomialForm -= currentP; // Kurangi untuk k[0]
                }
                BinomialForm += currentP; // Tambahkan untuk k[1]
            }
            
            cout << "Result P(" << k[0] << " <= X <= " << k[1]<< "): " << BinomialForm << endl;
            break;

        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
}

void HypergeometricFunction() {
    //Notasi N Populasi, K Elemen Sukses, n dan k cara
    int N, K, n, k;
    bool isBinomial;
    long double probability, HypergeomForm, mean, variance;

    do {
        cout << "1. Probability Mass Function" << endl;
        cout << "2. Mean and Variance Function" << endl;
        cout << "3. Cumulative Distributive Function" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option) {
        case 1:
            // Probability Formula
            cout << "Insert the amount of Sample (N) and Test Sample (n): ";
            cin >> N >> n;
            cout << "Insert the amount of Success (K) and Request Object to Find Probability (k): ";
            cin >> K >> k;
            
            HypergeomForm = (Combinationfact(K, k) * Combinationfact(N - K, n - k)) / Combinationfact(N, n);
            cout << "Probability Result: " << HypergeomForm << endl;
            break;

        case 2:
            // Mean and Variances Hypergeometric Distribution
            mean = 0.0, variance = 0.0;

            cout << endl << "Is this a Binomial Distribution? (1 for Yes, 0 for No): ";
            cin >> isBinomial;

            if (isBinomial) {
                cout << "Insert Mean and Variance: ";
                cin >> mean >> variance;

                // Calculate n and p for binomial distribution
                probability = (mean - variance) / mean;
                n = mean / probability;

                cout << "Binomial Parameters:" << endl;
                cout << "n (Number of Trials): " << n << endl;
                cout << "p (Probability of Success): " << probability << endl;
                BinomialFunction(probability);
            } else {
                cout << "Insert the amount of Sample (N) and Test Sample (n): ";
                cin >> N >> n;
                cout << "Insert the amount of Request Object to Find Probability (k): ";
                cin >> k;
                
                mean = n * (k / (double)N);
                variance = mean * (1 - (k / (double)N)) * ((double)(N - n) / (N - 1));
                
                cout << "Mean: " << mean << endl;
                cout << "Variance: " << variance << endl;
            }
            break;

        case 3:
            // Cumulative Distributive Function Hypergeometric Distribution
            cout << "Insert the amount of Sample (N) and Test Sample (n): ";
            cin >> N >> n;
            cout << "Insert the amount of Success (K) and Request Object to Find Probability (k): ";
            cin >> K >> k;

            HypergeomForm = 0.0;
            for (int i = 0; i <= k; i++){
                HypergeomForm += (Combinationfact(K, i) * Combinationfact(N - K, n - i)) / Combinationfact(N, n);
            }

            cout << "Result (X <= x): " << HypergeomForm << endl;
            cout << "Result (X >= x): " << 1 - HypergeomForm << endl;
            break;

        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
}

void PoissonFunction() {
    int n, k;
    bool isBinomial;
    long double probability, mean, PoissonForm;

    do {
        cout << "1. Probability Mass Function" << endl;
        cout << "2. Cumulative Distributive Function" << endl;
        cout << "3. Mean and Variance Function" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option){
        case 1:
            // Probability Mass Function
            mean = 0.0;
            cout << endl << "Is this a Binomial Distribution? (1 for Yes, 0 for No): ";
            cin >> isBinomial;

            if(isBinomial) {
                cout << "Insert the amount of Sample (n) and Success Probability: ";
                cin >> n >> probability;
                cout << "Insert the Request Object to Find that Probability: ";
                cin >> k;

                mean = n * probability;
            } else {
                cout << "Insert the Request Object to Find that Probability and Mean: ";
                cin >> k >> mean;
            }

            PoissonForm = (pow(mean, k) * pow(exp(1), -mean)) / factorial(k);;
            cout << "PMF Result: " << PoissonForm << endl;
            break;
        case 2:
            // Cumulative Distributive Function Poisson's Distribution
            mean = 0.0;
            cout << "Insert the Request Object to Find that Probability and Mean: ";
            cin >> k >> mean;

            PoissonForm = 0.0;
            for (int i = 0; i <= k; i++){
                PoissonForm += (pow(mean, i) * pow(exp(1), -mean)) / factorial(i);
            }

            cout << "Result (X <= x): " << PoissonForm << endl;
            cout << "Result (X >= x): " << 1 - PoissonForm << endl;
            break;
        case 3:
            // Mean and Variances Poisson's Distribution
            mean = 0.0;
            cout << "Insert the Request Object to Find that Probability and Mean: ";
            cin >> k >> mean;

            for (int i = 0; i <= k; i++){
                mean += i *(pow(mean, i) * pow(exp(1), -mean)) / factorial(i);;
            }

            // cout << mean
            break;
        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0)); 
}

void UniformFunction() {
    long double a, b, x1, x2, PDF, CDF, mean, variance;
    double diff;

    do {
        cout << "1. Probability Density Function" << endl;
        cout << "2. Cumulative Distributive Function" << endl;
        cout << "3. Mean and Variance Function" << endl;
        cout << "4. Finding interval [a, b]" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option) {
        case 1:
            // Probability Dense Function
            cout << "Input the interval A(Lower Limit) & B(Upper Limit): ";
            cin >> a >> b;
            cout << "Input the Request Interval X1 & X2: ";
            cin >> x1 >> x2;

            PDF = 0.0;
            if (a <= b && b >= x2 && x2 >= x1 && x1 >= a) {
                PDF = (1.0 / (b - a)) * (x2 - x1);
            } else {
                PDF = 0;
            }

            cout << "PDF Result: " << PDF << endl;
            break;
        case 2:
            // Cumulative Distributive Function
            double x;
            cout << "Input the interval A(Lower Limit) & B(Upper Limit): ";
            cin >> a >> b;
            cout << "Input the Request Interval X: ";
            cin >> x;

            CDF = 0.0;
            if (x <= a) CDF = 0;
            else if (a < x && x < b) {
                CDF = (x - a) / (b - a);
            } 
            else if (x >= b) CDF = 1;

            cout << "Cumulative Result: " << CDF << endl;
            break;
        case 3:
            // Mean and Variance Uniform Distribution
            cout << "Input the interval A(Lower Limit) & B(Upper Limit): ";
            cin >> a >> b;

            mean = 0.0, variance = 0.0;
            mean = (b + a) / 2;
            variance = pow((b - a), 2) / 12;

            cout << "Mean: " << mean << endl;
            cout << "Variance: " << variance << endl;
            break;
        case 4:
            // Finding Interval [a, b]
            cout << "Input the value of known Variable Mean and Variance: ";
            cin >> mean >> variance;

            // Calculate a and b
            diff = sqrt(12 * variance); // sqrt(12 * Variance)
            b = mean + (diff / 2.0);
            a = mean - (diff / 2.0);

            cout << "Calculated Interval: [" << a << ", " << b << "]" << endl;
            break;
        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
    
}

void ExponentFunction() {
    long double events, Rate_event, probability, mean, variance; // avg_time_in_event

    do {
        cout << "1. Probability Density Function" << endl;
        cout << "2. Cumulative Distributive Function" << endl;
        cout << "3. Mean and Variance" << endl;
        cout << "4. Finding missing value of λ Time" << endl;
        cout << "5. FInding missing value of Variance" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option) {
        case 1:
            // Probability Density Function
            cout << "Input the Time to represent the next event and Rate parameter: ";
            cin >> events >> Rate_event;

            if (events >= 0) {
                probability = Rate_event*pow(exp(1), -Rate_event*events);
            } else probability = 0;
            
            cout << "PDF Result: " << probability << endl;
            break;
        case 2:
            // Cumulative Distribution Function
            cout << "Input the Time to represent the next event and Rate parameter: ";
            cin >> events >> Rate_event;

            if (events >= 0) {
                probability = 1 - pow(exp(1), -Rate_event*events);
            } else probability = 0;
            
            cout << "CDF Result X >= x: " << 1 - probability << endl;
            cout << "CDF Result X <= x: " << probability << endl;
            break;
        case 3:
            // Mean and Variance
            cout << "Input the Rate Parameter: ";
            cin >> Rate_event;

            mean = 1 / Rate_event;
            variance = 1 / pow(Rate_event, 2);

            cout << "Mean: " << mean << endl;
            cout << "Variance: " << variance << endl;
            break;
        case 4:
            // Finding missing value of Event Time if known x and Probability
            cout << "Input the Probability and Rate Parameter: ";
            cin >> probability >> Rate_event;

            events = -(log10(probability) / log10(exp(1))) / Rate_event;

            cout << "Events(λ): " << events << endl;
            break;
        case 5:
            // Finding missing value of Variance if known equation in two Distribution
            cout << "Input the value of probability: ";
            cin >> probability;

            variance = 1.0 / sqrt(probability);

            cout << "Variance: " << variance << endl;
            break;
        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
}

void NormalDistFunction() {
    /* Friedrich Gauss (error theory), Abraham de Moivre (Central Limit Theorem)
       Pierre-Simon Laplace, Simeon Denis Poisson (Poisson's Distribution), Adolphe
       Quetelet the first Scientist of Sociolo  gy applying the Normal Distribution */
    long double x, Z_value, mean = 0.0, std_deviation = 0.0, variance, PDF, probability;

    do {
        cout << "1. Probability Density Function" << endl;
        cout << "2. Transform to Normal Standard-Z Score" << endl;
        cout << "3. Calculate Z-value Area Under the Curve" << endl;
        cout << "4. Calculate X by Z-Value" << endl;
        cout << "5. Mean and Variance" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option) {
        case 1:
            // Probability Density Function
            cout << "Input the value of Random Variable (x): ";
            cin >> x;
            cout << "Input the value of Mean and Std. Deviation: ";
            cin >> mean >> std_deviation;

            if (std_deviation <= 0) {
                cout << "Error: Standard deviation must be greater than 0!" << endl;
                break;
            }

            variance = pow(std_deviation, 2); // Update variance
            PDF = (1.0 / (std_deviation * sqrt(2 * M_PI))) * exp(-pow(x - mean, 2) / (2 * variance));

            cout << "PDF Result: " << PDF << endl;
            break;
        case 2:
            // Transforming Z into Normal Standard Score
            cout << "Input the value of Initial X: ";
            cin >> x;

            if (!mean && !std_deviation) {
                cout << "Input the value of Mean and Variance: ";
                cin >> mean >> std_deviation;
            }
            
            Z_value = (x - mean) / std_deviation;
            probability = 0.5 * (1 + erf(Z_value / sqrt(2)));

            cout << "Z-Value: " << Z_value << endl;
            cout << "Z-Value into Probability: " << Z_value << " into " << probability << endl;
            break;
        case 3:
            // Calculate Z-value Area if did/didn't known Z1 and Z2
            cout << "Input the value of Mean and Std. Deviation: ";
            cin >> mean >> std_deviation;
            cout << "Input the value of Initial X1 and X2: ";

            probability = 0.0;
            for (int i = 0; i < 2; i++) {
                cin >> x;
                Z_value = (x - mean) / std_deviation;
                cout << "Z" << i + 1 << ": " << Z_value << endl;

                // Probabilitas kumulatif untuk Z_value
                if (i == 0) {
                    probability += 0.5 * (1 + erf(Z_value / sqrt(2)));
                } else {
                    probability -= 0.5 * (1 + erf(Z_value / sqrt(2)));
                }
            }

            cout << "Result of P(a <= X <= b): " << probability << endl;
            break;
        case 4:
            // Calculate X if known probability, using Mean and Variance
            // Known as Inverse Gaussian Dist. or Wald Dist.
            cout << "Input the Probability (use dot if float): ";
            cin >> probability;
            cout << "Input the value of Mean and Std. Deviation: ";
            cin >> mean >> std_deviation;

            if (probability > 0 && probability < 1) {
                // Z_value = (1 / sqrt(2*M_PI*pow(probability, 3))) * exp(-(pow((probability - 1), 2) / 2*probability)); Inverse?
                Z_value = sqrt(2)*erf_inv(2 * probability - 1);
                cout << "Z-value: " << Z_value << endl;
                x = mean + (Z_value * std_deviation);
            } else cout << "Probability must be between 0 and 1 (exclusive)." << endl;

            cout << "Result P(X <= " << probability << "): " << x << endl;
            break;
        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
}

void MarkovChain() {
    int n, k, t, answer;
    long double stateProbability[10][10], statePi[10], sum[10] = {};
    
    do {
        cout << "1. Finding Matrixes Values" << endl;
        cout << "2. Calculate Transitive Probability Matrixes" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option) {
        case 1:
            // Finding Value of first state matrix
            cout << "Is there a Probability State(Pi) for next Periode? (1 for yes, other for no) ";
            cin >> answer;

            if (answer == 1) {
                cout << "Input the size of State Probability Matrixes: ";
                cin >> t;

                cout << "Input the Value of Pi(t): ";
                for (int i = 0; i < t; i++) {
                    cin >> statePi[i];
                }
                
                // Probability Matrixes given from Question
                cout << "Input the size of Transitive Probability (row, then column): ";
                cin >> k >> n;

                cout << "Input the Value: ";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        cin >> stateProbability[i][j];
                    }
                }
            } else {
                cout << "Input the size of Transitive Probability (row, then column): ";
                cin >> k >> n;

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        cin >> stateProbability[i][j];
                    }
                }
            }

            break;
        case 2:
            // Calculate State Probability (pi) with Transitive Probability
            if (answer == 1) {
                cout << "State (pi): ";
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        sum[i] += (statePi[j] * stateProbability[j][i]);
                    }
                    cout << sum[i] << " ";
                }

                cout << endl;
            } else {
                cout << "Input the Value to predict each moment: ";
                cin >> t;

                long double tempMatrix[10][10], resultMatrix[10][10];

                // Initialize resultMatrix as the identity matrix
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        resultMatrix[i][j] = (i == j) ? 1 : 0;
                    }
                }

                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        tempMatrix[i][j] = stateProbability[i][j];
                    }
                }
                
                for (int step = 0; step < t; step++) {
                    long double intermediateMatrix[10][10] = {};

                    // Perform matrix multiplication: resultMatrix * tempMatrix
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < k; j++) {
                            for (int l = 0; l < k; l++) {
                                intermediateMatrix[i][j] += resultMatrix[i][l] * tempMatrix[l][j];
                            }
                        }
                    }

                    // Update resultMatrix with the new values
                    for (int i = 0; i < n; i++) {
                        for (int j = 0; j < k; j++) {
                            resultMatrix[i][j] = intermediateMatrix[i][j];
                        }
                    }
                }

                // Print the result matrix P^t
                cout << "Probability Matrix after " << t << " step(s):" << endl;
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < k; j++) {
                        cout << resultMatrix[i][j] << " ";
                    }
                    cout << endl;
                }
            }
            
            break;
        default:
            break;
        }

        cout << "Done? ";
        cin.ignore();
        cin.getline(text, sizeof(text));
        system("cls");
    } while ((strcmp("No", text) == 0 || strcmp("no", text) == 0));
}

int main () {
    system("cls");
    long double prob;

    do {
        cout << "1. Bernouli's Distribution" << endl;
        cout << "2. Binomial Distribution" << endl;
        cout << "3. Hypergeometric Distribution" << endl;
        cout << "4. Poisson's Distribution" << endl;
        cout << "5. Uniform Distribution" << endl;
        cout << "6. Exponential Probability Distribution" << endl;
        cout << "7. Normal Distribution" << endl;
        cout << "8. Markov's Chain" << endl;
        cout << "Insert one of the option: ";
        cin >> option;
        system("cls");
        switch (option) {
        case 1:
            cout << "Insert the Probability Success (Default: Success): ";
            cin >> prob;
            BernouliFunction(prob, 1);
            break;
        case 2:
            cout << "Insert the Probability Success (Default: Success): ";
            cin >> prob;
            BinomialFunction(prob);
            break;
        case 3:
            HypergeometricFunction();
            break;
        case 4:
            PoissonFunction();
            break;
        case 5:
            UniformFunction();
            break;
        case 6:
            ExponentFunction();
            break;
        case 7:
            NormalDistFunction();
            break;
        case 8:
            MarkovChain();
            break;
        default:
            break;
        }
    } while (option >= 0 && option != 0);

    return 0;
}
