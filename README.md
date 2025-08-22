# orbital_propagator
Hardware Implementation using GNSS sensor and RPI pico of J2 orbital propagator for satellites in LEO.
Part of work I did as a member of the Student Satellite Programme, IITB

- Satellites for their functioning need to know where they are in space.
- GNSS sensors are typically used for this purpose.
- Constantly using GNSS sensors drains a lot of power of the satellite- a scarce resource- therefore we need an alternative.
- We periodically take GNSS readings
- In between successive readings, we use a mathematical model of the satellite kinematics to figure out its position.

The above has been implemented and validated on hardware(GNSS readings were not real but were emulated and parsed in the same format as they would normally be using a GNSS sensor)
