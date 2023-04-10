#include <iostream>
#include <fstream>
#include <glog/logging.h>
#include <eigen3/Eigen/Core>

const Eigen::Vector2d x0 = {10, 0};
const double u0 = 0;

double u = 0;
Eigen::Vector2d x = Eigen::Vector2d::Zero();

const double C = 18/M_PI;
const double ts = 0.2;
bool is_first = true;

Eigen::Vector2d calc_k(Eigen::Vector2d x, double u) {
    Eigen::Vector2d kx = Eigen::Vector2d::Zero();

    kx.x() = x.y();
    kx.y() = u - C*sin(x.x()/C);

    return kx;
}

Eigen::Vector2d rk4(Eigen::Vector2d x_k, double u_k, double dt) {
    Eigen::Vector2d x_k_1;
    Eigen::Vector2d k1, k2, k3, k4;
    k1 = calc_k(x_k, u_k);
    k2 = calc_k(x_k+0.5*dt*k1, u_k);
    k3 = calc_k(x_k+0.5*dt*k2, u_k);
    k4 = calc_k(x_k+dt*k3, u_k);

    x_k_1 = x_k + dt*(k1 + 2*k2 + 3*k3 + k4)/6;
    return x_k_1;
}

int main(int argc, char ** argv) {
    google::InitGoogleLogging(argv[0]);
    google::SetStderrLogging(google::INFO);

    std::ofstream ofs;
    ofs.open("rk4_data.csv", std::ios::out);
    if(!ofs.is_open()) {
        LOG(ERROR) << "open rk4_data.csv failed";
    }

    for(int i=0; i<50; i++) {
        if(is_first) {
            ofs << i*ts <<"," << x0(0) << "," << x0(1) << std::endl;
            x = rk4(x0, u0, ts);
            is_first = false;
        } else {
            x = rk4(x, u, ts);
        }
        ofs << i*ts << "," << x(0) << "," << x(1) << std::endl;
        LOG(INFO) << "x(0): " << x(0) << " x(1): " << x(1);
    }
    ofs.flush();
    ofs.close();
}
