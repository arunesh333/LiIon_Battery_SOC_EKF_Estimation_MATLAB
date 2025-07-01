# Li-ion Battery SOC Estimation using EKF (MATLAB)

A MATLAB project to estimate the **State of Charge (SOC)** of a lithium-ion battery using the **Extended Kalman Filter (EKF)** and a **first-order Thevenin model**. Built for **real-world EV battery data**, this simulation compares EKF results with traditional Coulomb Counting.

---

## What It Does

- Models Li-ion battery behavior with Thevenin RC circuit
- Fits OCV-SOC curve using 11th-order polynomial
- Applies EKF for SOC estimation using voltage, current, temperature inputs
- Compares accuracy against Coulomb Counting

---

## Results

- ✅ SOC error: ±1% (first 3 hrs), within ±2.5% overall  
- ✅ Voltage error: mostly within ±0.5 V  
- ✅ No drift over time; robust under noisy conditions

---

## Tools Used

- MATLAB
- Real EV battery dataset (`.mat`)
- Parameter interpolation using `scatteredInterpolant`

---

## Key Files

- `SOC_Estimation_by_EKF_method.m` — Main script  
- `BatteryModel.mat` & `SOC-OCV.mat` — Battery data  
- `Battery_datapoint.mat` — Input dataset

---

## Applications

- Battery Management Systems (BMS)
- EV SOC estimation
- Embedded battery control systems

---

##  Author

**Arunesh E** — Final Year EEE Student @ SASTRA  
Interested in EV systems & control modeling  
[LinkedIn](https://www.linkedin.com/in/arunesh33/)


---

