// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2018.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <cmath>
#include <limits>


double t = 0;
double tFinal = 0;
double tPlot = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double **x;

/**
 * Equivalent to x storing the velocities.
 */
double **v;

/**
 * One mass entry per molecule/particle.
 */
double *mass;

/**
 * Global time step size used.
 */
double timeStepSize = 0.000001;

/**
 * Maximum velocity of all particles.
 */
double maxV;

/**
 * Minimum distance between two elements.
 */
double minDx;


double **force  = new double *[NumberOfBodies];

/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char **argv) {
    NumberOfBodies = (argc - 2) / 7;

    x = new
            double *[NumberOfBodies];
    v = new
            double *[NumberOfBodies];
    mass = new
            double[NumberOfBodies];

    int readArgument = 1;

    tPlotDelta = std::stof(argv[readArgument]);
    readArgument++;
    tFinal = std::stof(argv[readArgument]);
    readArgument++;

    for (int i = 0; i < NumberOfBodies; i++) {
        x[i] = new
                double[3];
        v[i] = new
                double[3];

        x[i][0] = std::stof(argv[readArgument]);
        readArgument++;
        x[i][1] = std::stof(argv[readArgument]);
        readArgument++;
        x[i][2] = std::stof(argv[readArgument]);
        readArgument++;

        v[i][0] = std::stof(argv[readArgument]);
        readArgument++;
        v[i][1] = std::stof(argv[readArgument]);
        readArgument++;
        v[i][2] = std::stof(argv[readArgument]);
        readArgument++;

        mass[i] = std::stof(argv[readArgument]);
        readArgument++;

        if (mass[i] <= 0.0) {
            std::cerr << "invalid mass for body " << i << std::endl;
            exit(-2);
        }
    }

    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;

    if (tPlotDelta <= 0.0) {
        std::cout << "plotting switched off" << std::endl;
        tPlot = tFinal + 1.0;
    } else {
        std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
        tPlot = 0.0;
    }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
    videoFile.open("result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">"
              << std::endl
              << "<Collection>";
}


/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
    static int counter = -1;
    counter++;
    std::stringstream filename;
    filename << "result-" << counter << ".vtp";
    std::ofstream
            out(filename.str().c_str());
    out << "<VTKFile type=\"PolyData\" >" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
        << "  <Points>" << std::endl
        << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int i = 0; i < NumberOfBodies; i++) {
        out << x[i][0]
            << " "
            << x[i][1]
            << " "
            << x[i][2]
            << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>" << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>"
              << std::endl;
}


/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
    maxV = 0.0;
    minDx = std::numeric_limits<double>::max();

    int i, j;
    double xi, yi, zi, dx, dy, dz, F, fr2, fr6, fx, fy, fz, mt;

    double epsilon = 1.65e-21;
    double sigma = 3.4e-10;
    double sigma2 = sigma * sigma;

    // to avoid declaring the force for every run of UpdateBody it has been set on the class level
    // initialize the values for the forces 2D array
    for (i = 0; i < NumberOfBodies; i++) {
        force[i][0] = 0.0;
        force[i][1] = 0.0;
        force[i][2] = 0.0;
    }


    for (i = 0; i < NumberOfBodies; ++i) {
        xi = x[i][0];
        yi = x[i][1];
        zi = x[i][2];
        fx = 0.0;
        fy = 0.0;
        fz = 0.0;

        //reference for step2 http://phys.ubbcluj.ro/~tbeu/MD/C2_for.pdf
        // http://courses.cs.vt.edu/cs4414/S15/LECTURES/MolecularDynamics.pdf
        // http://phycomp.technion.ac.il/~talimu/md2.html
        // the last r is squared because we break the force down to x,y and z components
        for (j = 0; j < NumberOfBodies; j++) {
            if(i == j) continue;

            dx = xi - x[j][0];
            dy = yi - x[j][1];
            dz = zi - x[j][2];

            const double r2 = dx * dx + dy * dy + dz * dz;

            fr2 = sigma2 / r2;
            fr6 = fr2 * fr2 * fr2;

            /**
             * simplified Lennard-Jones
             * U = 4 * epsilon * ((sigma/r)^12 - (sigma/r)^6)
             * f = -dU/dr = 48 * epsilon * (sigma/r)^6 * ((sigma/r)^6 - 0.5) / r
             * fx = -dU/dx = -(x/r) * dU/dr
             *
             * finding the square root is an expensive operation so since all components of the force (x,y,z) are
             * divided by r, the force is divided by r2 as it is already computed from dx, dy and dz
             */

            F = 48.0 * epsilon * fr6 * (fr6 - 0.5) / r2;

            fx += dx * F;
            fy += dy * F;
            fz += dz * F;

            minDx = std::min(minDx, r2);
        }

        force[i][0] = fx;
        force[i][1] = fy;
        force[i][2] = fz;
    }

    minDx = std::sqrt(minDx);


    for (i = 0; i < NumberOfBodies; i++) {

        x[i][0] = x[i][0] + timeStepSize * v[i][0];
        x[i][1] = x[i][1] + timeStepSize * v[i][1];
        x[i][2] = x[i][2] + timeStepSize * v[i][2];

        mt = timeStepSize / mass[i];

        v[i][0] = v[i][0] + mt * force[i][0];
        v[i][1] = v[i][1] + mt * force[i][1];
        v[i][2] = v[i][2] + mt * force[i][2];

        const double V = std::sqrt(v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);

        maxV = std::max(maxV, V);
    }

    t += timeStepSize;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char **argv) {
    if (argc == 1) {
        std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time objects" << std::endl
                  << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting"
                  << std::endl
                  << "  final-time      simulated time (greater 0)" << std::endl
                  << std::endl
                  << "Examples:" << std::endl
                  << "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1"
                  << std::endl
                  << "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one"
                  << std::endl
                  << "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture"
                  << std::endl
                  << std::endl
                  << "In this naive code, only the first body moves" << std::endl;

        return -1;
    } else if ((argc - 3) % 7 != 0) {
        std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)"
                  << std::endl;
        return -2;
    }

    setUp(argc, argv);

    openParaviewVideoFile();

    int snapshotCounter = 0;
    if (t > tPlot) {
        printParaviewSnapshot();
        std::cout << "plotted initial setup" << std::endl;
        tPlot = tPlotDelta;
    }

    int timeStepCounter = 0;

    for (int i = 0; i < NumberOfBodies; ++i) {
        force[i] = new double[3];
    }

    while (t <= tFinal) {
        updateBody();
        timeStepCounter++;
        if (t >= tPlot) {
            printParaviewSnapshot();
            std::cout << "plot next snapshot"
                      << ",\t time step=" << timeStepCounter
                      << ",\t t=" << t
                      << ",\t dt=" << timeStepSize
                      << ",\t v_max=" << maxV
                      << ",\t dx_min=" << minDx
                      << std::endl;

            tPlot += tPlotDelta;
        }
    }

    closeParaviewVideoFile();

    return 0;
}
