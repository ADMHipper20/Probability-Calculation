#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <iomanip>

using namespace std;

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
    int a, b;
    long double PDF, mean, variance, probability;

    do {
        cout << "1. Probability Mass Function" << endl;
        cout << "2. Cumulative Distributive Function" << endl;
        cout << "3. Mean and Variance Function" << endl;
        cout << "Choose what to find: ";
        cin >> option;

        switch (option) {
        case 1:
            // Probability Dense Function
            break;
        case 2:
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
        default:
            break;
        }
    } while (option >= 0 && option != 0);

    return 0;
}