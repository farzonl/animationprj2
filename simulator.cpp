#include "simulator.h"
#include <fstream>
#include <iostream>
#define DEBUG 0
using namespace std;

double Equations::EARTH_GRAVITY = -9.81;

typedef void(*equationsFunctionPtr)(Particle &p,
                                  double timeStep);


void Equations::explicitEulerEquation(Particle &p,
                                      double timeStep) {
    
    Eigen::Vector3d acceleration = p.mAccumulatedForce / p.mMass;
    Eigen::Vector3d diplus1 = p.mPosition + p.mVelocity * timeStep;
    Eigen::Vector3d viplus1 = p.mVelocity + acceleration * timeStep;
    p.mPosition = diplus1;  // update position
    p.mVelocity = viplus1; // update veloicty
}

void Equations::midpointEquation(Particle &p,
                                 double timeStep) {

    Eigen::Vector3d acceleration = p.mAccumulatedForce / p.mMass;
    Eigen::Vector3d viplus1half = p.mVelocity + ((acceleration * timeStep) /2.0);
    Eigen::Vector3d diplus1 = p.mPosition + viplus1half * timeStep;
    Eigen::Vector3d viplus1 = p.mVelocity + acceleration * timeStep;
    p.mPosition = diplus1;  // update position
    p.mVelocity = viplus1; // update veloicty
}

void Equations::implicitEulerEquation(Particle &p,
                                  double timeStep) {
    Eigen::Vector3d acceleration = p.mAccumulatedForce / p.mMass;
    Eigen::Vector3d viplus1 = p.mVelocity + acceleration * timeStep;
    Eigen::Vector3d diplus1 = p.mPosition + viplus1 * timeStep;
    
    p.mPosition = diplus1;  // update position
    p.mVelocity = viplus1; // update veloicty
}

void ParticleSystem::initW(std::vector<Particle> &particles) {
    size_t matSize = 3 * particles.size();
    Eigen::MatrixXf W(matSize, matSize);
    W << Eigen::MatrixXf::Zero(matSize, matSize);

    for(int i = 0, j =0; i < particles.size(); i++, j+=3) {
        double inverseMass = 1.0 / particles[i].mMass;
        W(j,j) = inverseMass;
        W(j+1,j+1) = inverseMass;
        W(j+2,j+2) = inverseMass;
    }
    this->W = W;
}

Eigen::MatrixXf ParticleSystem::genJMat(std::vector<Particle> &particles, DERV_LEVEL dl) {
    size_t matSize = 3 * particles.size();
    Eigen::MatrixXf J(particles.size(), matSize);
    J << Eigen::MatrixXf::Zero(particles.size(),matSize);

    Eigen::Vector3d currVector;
    Eigen::Vector3d nextVector;
    // 2 constriants increment by 2 
    for(int i = 0; i < particles.size(); i+=2) {
        switch(dl) {
            case NO_DOT:
                currVector = particles[i].mPosition;
                nextVector = particles[i + 1].mPosition;
                break;
            case DOT:
                currVector = particles[i].mVelocity;
                nextVector = particles[i + 1].mVelocity;
                break;
        }
        // (1) loop constraint
        J(i, 0) = currVector[0];
        J(i, 1) = currVector[1];
        J(i, 2) = currVector[2];

        // (2) fixed distant contraint
        J(i+1, 0) = currVector[0] - nextVector[0];
        J(i+1, 1) = currVector[1] - nextVector[1];
        J(i+1, 2) = currVector[2] - nextVector[2];

        J(i+1, 3) = nextVector[0] - currVector[0];
        J(i+1, 4) = nextVector[1] - currVector[1];
        J(i+1, 5) = nextVector[2] - currVector[2];
    }

    //std::cout << "calc J: \n" << J << std::endl;
    return J;
}
void ParticleSystem::calcQ(std::vector<Particle> &particles) {
    size_t matSize = 3 * particles.size();
    Eigen::MatrixXf Q(matSize, 1);
    Q << Eigen::MatrixXf::Zero(matSize,1);
    for(int i = 0, j = 0; i < particles.size(); i++,j+=3) {
        Q(j,0) = particles[i].mAccumulatedForce[0];
        Q(j+1,0) = particles[i].mAccumulatedForce[1];
        Q(j+2,0) = particles[i].mAccumulatedForce[2];
    }
    this->Q = Q;
}

void ParticleSystem::calcq_dot(std::vector<Particle> &particles) {
    size_t matSize = 3 * particles.size();
    Eigen::MatrixXf q_dot(matSize, 1);
    q_dot << Eigen::MatrixXf::Zero(matSize,1);
    for(int i = 0, j = 0; i < particles.size(); i++,j+=3) {
        q_dot(j, 0) = particles[i].mVelocity[0];
        q_dot(j+1,0) = particles[i].mVelocity[1];
        q_dot(j+2,0) = particles[i].mVelocity[2];
    }
    this->q_dot = q_dot;
}


void ParticleSystem::calcLambda(std::vector<Particle> &particles) {
    // JWJ^T lambda = -J_dot * q_dot - JWQ
    // lambda = (-J_dot * qdot - JWQ)/JWJ^T
    this->J = ParticleSystem::genJMat(particles);
    //std::cout << "calc J: \n" << this->J << std::endl; 

    this->J_dot = ParticleSystem::genJMat(particles, DOT);
    //std::cout << "calc J_dot: \n" << J_dot << std::endl;

    calcQ(particles);
    calcq_dot(particles);

    //std::cout << "calc Q: \n" << Q << std::endl; 
    //std::cout << "calc q_dot: \n" << q_dot << std::endl;

    ///-J_dot * q_dot - JWQ
    Eigen::MatrixXf  JWQ = this->J * this->W * this->Q;
    //std::cout << "calc JWQ: \n" << JWQ << std::endl;
    Eigen::MatrixXf rightHandSide = (-1 * this->J_dot * this->q_dot - JWQ);
    //std::cout << "calc rightHandSide: \n" << rightHandSide << std::endl;
    // 1/JWJ^T = (JWJ^T)^-1
    Eigen::MatrixXf leftHandInverse = (this->J * this->W  * this->J.transpose()).inverse();
    //std::cout << "calc leftHandInverse: \n" << leftHandInverse << std::endl;
    this->lambda = leftHandInverse * rightHandSide ;
    //std::cout << "calc lambda: \n" << lambda << std::endl; 
    
}

void ParticleSystem::calcQhat(std::vector<Particle> &particles) {
    calcLambda(particles);
    //Once the linear system has been solved, the vector λ is 
    // multiplied by J^T to produce the global constraint force vector Q_hat
    this->Q_hat = J.transpose() * this->lambda;
    //std::cout << "calc Q_hat: \n" << Q_hat << std::endl;
}

void ParticleSystem::diagnostics() {
#if DEBUG
    std::cout << "calc W: \n" << W << std::endl;
    std::cout << "calc J: \n" << J << std::endl;
    std::cout << "calc J_dot: \n" << J_dot<< std::endl;
    std::cout << "calc Q: \n" << Q << std::endl;
    std::cout << "calc q_dot: \n" << q_dot<< std::endl;
    std::cout << "calc lambda: \n" << lambda << std::endl;
    std::cout << "calc Q_hat: \n" << Q_hat << std::endl;
#endif
}

void ParticleSystem::apply(std::vector<Particle> &particles, float timeStep,IntegraionMethod im) {
    equationsFunctionPtr fPtr = NULL;
    switch(im) {
        case explicitEuler:
            fPtr = &Equations::explicitEulerEquation;
            break;
        case midpoint:
            fPtr = &Equations::midpointEquation;
            break;
        case implicitEuiler:
            fPtr = &Equations::implicitEulerEquation;
            break;
    }
    
    calcQhat(particles);
    diagnostics();
    //std::cout << "P1 mAccumulatedForce before Q_hat: \n" << particles[0].mAccumulatedForce << std::endl;
    //std::cout << "P2 mAccumulatedForce before Q_hat: \n" << particles[2].mAccumulatedForce << std::endl;
    
    for(int i = 0, j = 0; i < particles.size(); i++, j+=3) {
        particles[i].mAccumulatedForce[0] += this->Q_hat(j, 0);
        particles[i].mAccumulatedForce[1] += this->Q_hat(j+1, 0);
        particles[i].mAccumulatedForce[2] += this->Q_hat(j+2, 0);
    }

     //std::cout << "P1 mAccumulatedForce after Q_hat: \n" << particles[0].mAccumulatedForce << std::endl;
     //std::cout << "P2 mAccumulatedForce after Q_hat: \n" << particles[2].mAccumulatedForce << std::endl;

    for(int i = 0; i < particles.size(); i++) {
        fPtr(particles[i], timeStep);
    }
}

Simulator::Simulator() : ps() {
    // initialize the particles
    mParticles.resize(2);
    reset();
    std::cout <<"####### START CONSTRUCTOR ###########" << std::endl;
    std::cout << "P1 mAccumulatedForce before: \n" << mParticles[0].mAccumulatedForce << std::endl;
    std::cout << "P2 mAccumulatedForce before: \n" << mParticles[2].mAccumulatedForce << std::endl;
    std::cout << "Particle velocity 1 before: \n" << mParticles[0].mVelocity << std::endl;
    std::cout << "Particle velocity 2 before: \n"  << mParticles[1].mVelocity << std::endl;
    std::cout << "Particle pos 1 before: \n" << mParticles[0].mPosition << std::endl;
    std::cout << "Particle pos 2 before: \n"  << mParticles[1].mPosition << std::endl;
    std::cout <<"####### END CONSTRUCTOR ###########" << std::endl;

    ps.initW(mParticles);
    //std::cout << "calc W: \n" << ps.W << std::endl;
    mTimeStep = 0.0003;
}

int Simulator::getNumParticles() {
    return mParticles.size();
}

Particle* Simulator::getParticle(int index) {
    return &mParticles[index];
}

double Simulator::getTimeStep() {
    return mTimeStep;
}

void Simulator::reset() {
    mParticles[0].mPosition[1] = 0.0;
    mParticles[0].mPosition[0] = 0.2;
    mParticles[1].mPosition[0] = 0.2;
    mParticles[1].mPosition[1] = -0.1;
    
    for (int i = 0; i < getNumParticles(); i++) {
        mParticles[i].mVelocity.setZero();
        mParticles[i].mAccumulatedForce.setZero();
    }
    
}

void Simulator::simulate() {
    //std::cout << "P1 mAccumulatedForce before: \n" << mParticles[0].mAccumulatedForce << std::endl;
    //std::cout << "P2 mAccumulatedForce before: \n" << mParticles[2].mAccumulatedForce << std::endl;
    //std::cout << "Particle velocity 1 before: \n" << mParticles[0].mVelocity << std::endl;
    //std::cout << "Particle velocity 2 before: \n"  << mParticles[1].mVelocity << std::endl;
    //std::cout << "Particle pos 1 before: \n" << mParticles[0].mPosition << std::endl;
    //std::cout << "Particle pos 2 before: \n"  << mParticles[1].mPosition << std::endl;
    for (int i = 0; i < mParticles.size(); i++) {
        mParticles[i].mAccumulatedForce[1] += Equations::EARTH_GRAVITY * mParticles[i].mMass;
    }
    ps.apply(mParticles,mTimeStep);

    //std::cout << "P1 mAccumulatedForce after: \n" << mParticles[0].mAccumulatedForce << std::endl;
    //std::cout << "P2 mAccumulatedForce after: \n" << mParticles[2].mAccumulatedForce << std::endl;
    //std::cout << "Particle velocity 1 after: \n" << mParticles[0].mVelocity << std::endl;
    //std::cout << "Particle velocity 2 after: \n"  << mParticles[1].mVelocity << std::endl;
    //std::cout << "Particle pos 1 after: \n" << mParticles[0].mPosition << std::endl;
    //std::cout << "Particle pos 2 after: \n"  << mParticles[1].mPosition << std::endl;
    
    for (int i = 0; i < mParticles.size(); i++) {
        mParticles[i].mAccumulatedForce.setZero();
    }
}







