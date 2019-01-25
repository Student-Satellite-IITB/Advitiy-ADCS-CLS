import numpy as np

#------------Initial conditions
v_q0_BO = np.array([0.,0.,0.,1.])	#unit quaternion initial condition
v_w0_BOB = np.array([0.,0.,0.])

MODEL_STEP=0.1
CONTROL_STEP = 2.0	#control cycle time period in second

# boolean determining whether satellite in SSO or PO
# 1 means SSO
# 0 means PO
orbitbool = 1

# boolean determining whether disturbance will act or not
# 1 means disturbance present
# 0 means disturbance absent
distbool = 0

# boolean determining whether sensor modelling will be used or not
# 1 means sensor-modelling present
# 0 means sensor modelling absent
sensbool = 0

# boolean determining whether estimator will be used or not
# 1 means estimator used
# 0 means estimator not used
estbool = 0

# constant determining which controller will act or not
# 1 means detumbing_controller present
# 0 means controller absent
contcons = 0


# boolean determining whether actuator modelling will be used or not
# 1 means actuator-modelling present
# 0 means actuator modelling absent
actbool = 0