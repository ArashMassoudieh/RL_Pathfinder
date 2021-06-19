#ifndef REINFORCEMENTLEARNER_H
#define REINFORCEMENTLEARNER_H

struct parameters
{
    double Gamma = 0.9;
    double epsilon = 0.1;
    double UltimateReward = 100000;
};

class ReinforcementLearner
{
    public:
        ReinforcementLearner();
        virtual ~ReinforcementLearner();
        ReinforcementLearner(const ReinforcementLearner& other);
        ReinforcementLearner& operator=(const ReinforcementLearner& other);
        parameters Parameters;
    protected:

    private:
};

#endif // REINFORCEMENTLEARNER_H
