#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <string>
#include "pico/stdlib.h"
#include "hardware/uart.h"
#define WGS84_A 6378137.0      // Semi-major axis in meters
#define WGS84_E2 0.00669437999014  // Square of eccentricity
#define DEG_TO_RAD(deg) ((deg) * (3.14159 / 180.0))
#define MU 398600.4418 // Standard gravitational parameter for Earth, km^3/s^2
#define J2000_EPOCH 2451545.0 // Julian date for J2000 epoch
#define SECONDS_IN_A_DAY 86400.0
#define UART_ID uart1
#define BAUD_RATE 115200
#define UART_TX_PIN 4
#define UART_RX_PIN 5
/*void writeStateToCSV(const std::string &filename, double ts, const double state_vector[6]) {
    std::ofstream file(filename, std::ios::app);
    file << ts;
    file << ",";
    for (int i = 0; i < 6; ++i) {
        file << state_vector[i];
        if (i < 5) {
            file << ",";  // Add a comma between elements, except at the end of the row
        }
    }
    file << "\n";  // Newline at the end of each row

    file.close();
}*/
double newtonraphson(double m, double e) {
	double ecan;
	if (m < 3.14159) {
		ecan = m + e / 2;
	} else {
		ecan = m - e / 2;
	}
	double f = ecan - e * sin(ecan) - m;
	double fdash = 1 - e * cos(ecan);
	double ratio = f / fdash;
	double tolerance = 1e-8;
	while (fabs(ratio) > tolerance) {
		ecan = ecan - ratio;
		f = ecan - e * sin(ecan) - m;
		fdash = 1 - e * cos(ecan);
		ratio = f / fdash;
	}
	return ecan;
}
void oe_time(double oe[], double deltat) {

	double e = oe[1];
	double k = sqrt((1 + e) / (1 - e));
	double theta0 = oe[5];

	double eccan0 = 2 * atan((1 / k) * tan(theta0 / 2));
	double M0 = eccan0 - e * sin(eccan0);

	double a = oe[0];
	double timep = 2 * 3.14159 * pow(a, 1.5) / sqrt(398600);
	double t0 = timep * M0 / (2 * 3.14159);

	double M = 2 * 3.14159 * (t0 + deltat) / timep;

	double eccan = newtonraphson(M, e);

	double theta = 2 * atan(k * tan(eccan / 2));
        oe[5]=theta;

}
double convertLatToDecimalDegrees(char *coordinate, char *direction) {
    double decimalDegrees = 0.0;
    double degrees, minutes;
    int len = strlen(coordinate);

    // Degrees are the first two digits
    char degreesStr[3];
    strncpy(degreesStr, coordinate, 2);
    degreesStr[2] = '\0';

    degrees = atof(degreesStr);

    // Minutes are the rest
    minutes = atof(coordinate + 2) / 60.0;

    decimalDegrees = degrees + minutes;

    if (direction[0] == 'S' || direction[0] == 'W') {
        decimalDegrees *= -1;
    }

    return decimalDegrees;
}
double convertLongToDecimalDegrees(char *coordinate, char *direction) {
    double decimalDegrees = 0.0;
    double degrees, minutes;
    int len = strlen(coordinate);

    // Degrees are the first two digits
    char degreesStr[4];
    strncpy(degreesStr, coordinate, 3);
    degreesStr[3] = '\0';

    degrees = atof(degreesStr);

    // Minutes are the rest
    minutes = atof(coordinate + 3) / 60.0;

    decimalDegrees = degrees + minutes;

    if (direction[0] == 'S' || direction[0] == 'W') {
        decimalDegrees *= -1;
    }

    return decimalDegrees;
}
void parseGPGGA(char *gpgga, double *latitude, double *longitude, double *altitude, double *timeInSeconds, double *hs, double *ms, double *ss) {
    char *token;
    int fieldIndex = 0;
    token = strtok(gpgga, ",");

    while (token != NULL) {
        switch (fieldIndex) {

            case 1: // Time
                {
                    int hours = (token[0] - '0') * 10 + (token[1] - '0');
                    int minutes = (token[2] - '0') * 10 + (token[3] - '0');
                    double seconds = (token[4] - '0') * 10 + (token[5] - '0') + (token[7] - '0')*0.1 + (token[8] - '0')*0.01;
                    *timeInSeconds = hours * 3600.0 + minutes * 60.0 + seconds;
                    *hs = hours*1.0;
                    *ms = minutes*1.0;
                    *ss = seconds*1.0;
                }
                break;
            case 2: // Latitude
                *latitude = convertLatToDecimalDegrees(token, strtok(NULL, ","));
                fieldIndex++; // Skip the direction field
                break;
            case 4: // Longitude
                *longitude = convertLongToDecimalDegrees(token, strtok(NULL, ","));
                fieldIndex++; // Skip the direction field
                break;
            case 9: // Altitude
                *altitude = atof(token);
                break;
            default:
                break;
        }
        token = strtok(NULL, ",");
        fieldIndex++;
    }
}
double utc_to_julian(int year, int month, int day, int hour, int minute, int second) {
    if (month <= 2) {
        year--;
        month += 12;
    }
    int A = year / 100;
    int B = 2 - A + A / 4;
    double julian_date = (int)(365.25 * (year + 4716)) + (int)(30.6001 * (month + 1)) + day + B - 1524.5;
    julian_date += (hour + minute / 60.0 + second / 3600.0) / 24.0; 
    return julian_date;
}
double calculate_gmst(double jd) {
    double t = (jd - J2000_EPOCH) / 36525.0; 
    double gmst = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t - 6.2e-6 * t * t * t; 
    gmst = fmod(gmst, SECONDS_IN_A_DAY); 
    if (gmst < 0) gmst += SECONDS_IN_A_DAY; 
    return gmst * (2 * 3.14159 / SECONDS_IN_A_DAY);
}
void ECEFtoECI(double ecef_x, double ecef_y, double ecef_z, 
                 int year, int month, int day, int hour, int minute, int second, 
                 double *eci_x, double *eci_y, double *eci_z) {
    double jd = utc_to_julian(year, month, day, hour, minute, second);
    double theta = calculate_gmst(jd);

    *eci_y = -ecef_x * cos(theta) - ecef_y * sin(theta);
    *eci_x = -ecef_x * sin(theta) + ecef_y * cos(theta);
    *eci_z = ecef_z; 
}
void coe_from_sv(const double R[3], const double V[3], double mu, double coe[6]) {
    const double eps = 1.e-10;
    
    double r = std::sqrt(R[0] * R[0] + R[1] * R[1] + R[2] * R[2]);
    double v = std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
    double vr = (R[0] * V[0] + R[1] * V[1] + R[2] * V[2]) / r;

    // Angular momentum vector and its magnitude
    double H[3] = {
        R[1] * V[2] - R[2] * V[1],
        R[2] * V[0] - R[0] * V[2],
        R[0] * V[1] - R[1] * V[0]
    };
    double h = std::sqrt(H[0] * H[0] + H[1] * H[1] + H[2] * H[2]);

    // Inclination (i)
    coe[2] = std::acos(H[2] / h);

    // Node line vector and its magnitude
    double N[3] = {-H[1], H[0], 0};
    double n = std::sqrt(N[0] * N[0] + N[1] * N[1]);

    // Right ascension of the ascending node (RA)
    coe[4] = (n != 0) ? std::acos(N[0] / n) : 0;
    if (n != 0 && N[1] < 0) coe[4] = 2 * 3.14159 - coe[4];

    // Eccentricity vector and its magnitude (e)
    double E[3] = {
        (v * v - mu / r) * R[0] - r * vr * V[0],
        (v * v - mu / r) * R[1] - r * vr * V[1],
        (v * v - mu / r) * R[2] - r * vr * V[2]
    };
    for (int i = 0; i < 3; ++i) E[i] /= mu;
    coe[1] = std::sqrt(E[0] * E[0] + E[1] * E[1] + E[2] * E[2]);

    // Argument of perigee (w)
    coe[3] = 0;
    if (n != 0 && coe[1] > eps) {
        coe[3] = std::acos((N[0] * E[0] + N[1] * E[1]) / (n * coe[1]));
        if (E[2] < 0) coe[3] = 2 * 3.14159 - coe[3];
    }

    // True anomaly (theta)
    coe[5] = 0;
    if (coe[1] > eps) {
        coe[5] = std::acos((E[0] * R[0] + E[1] * R[1] + E[2] * R[2]) / (coe[1] * r));
        if (vr < 0) coe[5] = 2 * 3.14159 - coe[5];
    } else {
        double cp = N[1] * R[2] - N[2] * R[1];
        coe[5] = (cp >= 0) ? std::acos((N[0] * R[0] + N[1] * R[1]) / (n * r)) : 2 * 3.14159 - std::acos((N[0] * R[0] + N[1] * R[1]) / (n * r));
    }

    // Semi-major axis (a)
    coe[0] = h * h / mu / (1 - coe[1] * coe[1]);
}
void LLAtoECEF(double lat, double lon, double alt, double *x, double *y, double *z) {
    double N; // Prime vertical radius of curvature
    double sinLat = sin(DEG_TO_RAD(lat));
    double cosLat = cos(DEG_TO_RAD(lat));
    double sinLon = sin(DEG_TO_RAD(lon));
    double cosLon = cos(DEG_TO_RAD(lon));

    N = WGS84_A / sqrt(1 - WGS84_E2 * sinLat * sinLat);

    *x = (N + alt) * cosLat * cosLon;
    *y = (N + alt) * cosLat * sinLon;
    *z = ((1 - WGS84_E2) * N + alt) * sinLat;
}
void matrix_multiply(const double A[3][3], const double B[3][3], double C[3][3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < 3; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Compute the rotation matrix
void compute_rotation_matrix(double Omega, double inc, double omega, double R[3][3]) {
    // Rotation matrix components
    double R_Omega[3][3] = {
        {cos(Omega), -sin(Omega), 0},
        {sin(Omega), cos(Omega), 0},
        {0, 0, 1}
    };

    double R_inc[3][3] = {
        {1, 0, 0},
        {0, cos(inc), -sin(inc)},
        {0, sin(inc), cos(inc)}
    };

    double R_omega[3][3] = {
        {cos(omega), -sin(omega), 0},
        {sin(omega), cos(omega), 0},
        {0, 0, 1}
    };

    double temp[3][3];
    matrix_multiply(R_Omega, R_inc, temp); // temp = R_Omega * R_inc
    matrix_multiply(temp, R_omega, R);     // R = temp * R_omega
}

// Main function to compute state vector from orbital elements
void sv_from_coe(const double coe[6], double state_vector[6]) {
    double a = coe[0];
    double e = coe[1];
    double inc = coe[2];
    double omega = coe[3];
    double Omega = coe[4];
    double theta = coe[5];

    // Semi-latus rectum
    double p = a * (1 - e * e);

    // Radius
    double r = p / (1 + e * cos(theta));

    // Position in perifocal (PQW) frame
    double r_pqw[3] = {r * cos(theta), r * sin(theta), 0};

    // Velocity in perifocal (PQW) frame
    double v_pqw[3] = {-sqrt(398600 / p) * sin(theta),
                       sqrt(398600 / p) * (e + cos(theta)),
                       0};

    // Rotation matrix from PQW to ECI
    double R[3][3];
    compute_rotation_matrix(Omega, inc, omega, R);

    // Convert position and velocity to ECI frame and store in state_vector
    for (int j = 0; j < 3; ++j) {
        state_vector[j] = R[j][0] * r_pqw[0] + R[j][1] * r_pqw[1] + R[j][2] * r_pqw[2];      // Position part
        state_vector[j + 3] = R[j][0] * v_pqw[0] + R[j][1] * v_pqw[1] + R[j][2] * v_pqw[2];  // Velocity part
    }
}
#define MAX_STRINGS 5         // Number of strings per buffer
#define STRING_LENGTH 100     // Length of each string (including null terminator)

char gpggaBuffer1[MAX_STRINGS][STRING_LENGTH];
char gpggaBuffer2[MAX_STRINGS][STRING_LENGTH];
char tempBuffer[STRING_LENGTH];  // Temporary buffer for incoming string
double latitude1, longitude1, altitude1, timeSeconds1, hs1, ms1, ss1;
double x1, y01, z1, x_eci1, y_eci1, z_eci1;
double latitude2, longitude2, altitude2, timeSeconds2, hs2, ms2, ss2;
double x2, y2, z2, x_eci2, y_eci2, z_eci2;
double r_roll_5[3], v_roll_5[3];
double r[3];
double v[3];
void store_in_buffer(char* source, int index, char targetBuffer[MAX_STRINGS][STRING_LENGTH]) {
    for (int i = 0; i < STRING_LENGTH - 1 && source[i] != '\0'; i++) {
        targetBuffer[index][i] = source[i];
    }
    targetBuffer[index][STRING_LENGTH - 1] = '\0';  // Null-terminate the string
}
void parse_strings(){
    for(int roll=0; roll<5; roll++){
        parseGPGGA(gpggaBuffer1[roll], &latitude1, &longitude1, &altitude1, &timeSeconds1, &hs1, &ms1, &ss1);
        parseGPGGA(gpggaBuffer2[roll], &latitude2, &longitude2, &altitude2, &timeSeconds2, &hs2, &ms2, &ss2);
        LLAtoECEF(latitude1, longitude1, altitude1, &x1, &y01, &z1);
        ECEFtoECI(x1, y01, z1, 2024, 8, 6, hs1, ms1, ss1, &x_eci1, &y_eci1, &z_eci1);//date from RTC
        LLAtoECEF(latitude2, longitude2, altitude2, &x2, &y2, &z2);
        ECEFtoECI(x2, y2, z2, 2024, 8, 6, hs2, ms2, ss2, &x_eci2, &y_eci2, &z_eci2);//date from RTC
        double v1, v2, v3;
        v1=(x_eci2-x_eci1)/(timeSeconds2-timeSeconds1);
        v2=(y_eci2-y_eci1)/(timeSeconds2-timeSeconds1);
        v3=(z_eci2-z_eci1)/(timeSeconds2-timeSeconds1);
        r[0]=x_eci2/1000;
        r[1]=y_eci2/1000;
        r[2]=z_eci2/1000;
        v[0]=v1/1000;
        v[1]=v2/1000;
        v[2]=v3/1000;
        r_roll_5[0]=(r_roll_5[0]*roll+r[0])/(roll+1);
        r_roll_5[1]=(r_roll_5[1]*roll+r[1])/(roll+1);
        r_roll_5[2]=(r_roll_5[2]*roll+r[2])/(roll+1);
        v_roll_5[0]=(v_roll_5[0]*roll+v[0])/(roll+1);
        v_roll_5[1]=(v_roll_5[1]*roll+v[1])/(roll+1);
        v_roll_5[2]=(v_roll_5[2]*roll+v[2])/(roll+1);
    }
    double timeSeconds=0.5*(timeSeconds1+timeSeconds2);
    double coe[6];
    coe_from_sv(r_roll_5, v_roll_5, 398600, coe);
    printf("Initial OE and State\n");
    printf("%f\n", coe[0]);
    printf("%f\n", coe[1]);
    printf("%f\n", coe[2]);
    printf("%f\n", coe[3]);
    printf("%f\n", coe[4]);
    printf("%f\n", coe[5]);
    printf("\n");
    printf("%f\n", r_roll_5[0]);    
    printf("%f\n", r_roll_5[1]);
    printf("%f\n", r_roll_5[2]);
    printf("%f\n", v_roll_5[0]);    
    printf("%f\n", v_roll_5[1]);
    printf("%f\n", v_roll_5[2]);
    /*std::string filename = "state_vectors.csv";
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    file << "t,X,Y,Z,Vx,Vy,Vz\n"; 
    file.close();*/
    for(int i=0; i<500; i++){
        oe_time(coe, 1);
        double sv[6];
        sv_from_coe(coe, sv);
        printf("OE and state at %f seconds\n", timeSeconds+i);
        printf("%f\n", coe[0]);
        printf("%f\n", coe[1]);
        printf("%f\n", coe[2]);
        printf("%f\n", coe[3]);
        printf("%f\n", coe[4]);
        printf("%f\n", coe[5]);
        printf("\n");
        printf("%f\n", sv[0]);
        printf("%f\n", sv[1]);
        printf("%f\n", sv[2]);
        printf("%f\n", sv[3]);
        printf("%f\n", sv[4]);
        printf("%f\n", sv[5]);
        //writeStateToCSV(filename, timeSeconds+i, sv);
    }
}
int bufferIndex = 0;
char gpggaBuffer[100];
int counter=0;


int main() {
    stdio_init_all();
    uart_init(UART_ID, BAUD_RATE);
    gpio_set_function(UART_TX_PIN, GPIO_FUNC_UART);
    gpio_set_function(UART_RX_PIN, GPIO_FUNC_UART);

    printf("UART initialized. Starting to read GPGGA strings...\n");

    while (true) {
        if (uart_is_readable(UART_ID)) {
            char incomingByte = uart_getc(UART_ID);

            if (incomingByte == '$') {
                // Start of a new GPGGA message
                bufferIndex = 0;
            }

            if (bufferIndex < STRING_LENGTH - 1) {
                tempBuffer[bufferIndex++] = incomingByte;
            }

            if (incomingByte == '\n') {
                // End of message; null-terminate and store
                tempBuffer[bufferIndex - 1] = '\0';  // Remove newline and null-terminate

                if (counter % 2 == 0) {
                    store_in_buffer(tempBuffer, counter / 2, gpggaBuffer1);
                    printf("Stored in gpggaBuffer1[%d]: %s\n", counter / 2, gpggaBuffer1[counter / 2]);
                } else {
                    store_in_buffer(tempBuffer, (counter - 1) / 2, gpggaBuffer2);
                    printf("Stored in gpggaBuffer2[%d]: %s\n", (counter - 1) / 2, gpggaBuffer2[(counter - 1) / 2]);
                }

                counter++;
                if (counter >= MAX_STRINGS * 2) {
                    // Run parse_strings() after filling both buffers
                    parse_strings();
                    break;  // Exit the loop after parsing
                }

                bufferIndex = 0;  // Reset temp buffer index for next message
            }
        }
    }

    printf("Process completed.\n");
    return 0;
}