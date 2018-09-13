#ifndef SIMULATOR_H
#define SIMULATOR_H

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT

#include <Eigen/Dense>
#include <vector>
#include <string>

#include "particle.h"

class Equations {
public:
    static double EARTH_GRAVITY;

    static void explicitEulerEquation(Particle &p,
                                  double timeStep);
    
    static void midpointEquation(Particle &p,
                                      double timeStep);
    
    static void implicitEulerEquation(Particle &p,
                                      double timeStep);
private:
    Equations() {}
       
};

enum IntegraionMethod {
    explicitEuler,
    midpoint,
    implicitEuiler
};

enum DERV_LEVEL{
    NO_DOT, // pos
    DOT    // velocity
    //,DOT_DOT // acceleration
};

class ParticleSystem{
    public:
        void initW(std::vector<Particle> &particles);
        void apply(std::vector<Particle> &particles, float timeStep,IntegraionMethod im=midpoint);
        Eigen::MatrixXf W;
        Eigen::MatrixXf J;
        Eigen::MatrixXf J_dot;
        Eigen::MatrixXf Q; 
        Eigen::MatrixXf q_dot;
        Eigen::MatrixXf lambda;
        Eigen::MatrixXf Q_hat;
    private:
        static Eigen::MatrixXf genJMat(std::vector<Particle> &particles, DERV_LEVEL dl=NO_DOT);
        void calcQhat(std::vector<Particle> &particles);
        void calcLambda(std::vector<Particle> &particles);
        void calcQ(std::vector<Particle> &particles);
        void calcq_dot(std::vector<Particle> &particles);
        void diagnostics();
};

// class containing objects to be simulated
class Simulator {
public:
    Simulator();
        
    void simulate();
    
    int getNumParticles();
    
    Particle* getParticle(int);
    
    double getTimeStep();
    
    void reset();
private:
    double mTimeStep;       // time step
    std::vector<Particle> mParticles;
    ParticleSystem ps;
    
};

#endif  // SIMULATOR_H
