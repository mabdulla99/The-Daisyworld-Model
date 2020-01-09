############################## Defining functions for Daisyworld with plagues

def plague(t): 
  return gamma + 0.6*abs(np.sin(np.pi*t/100))

def alphaw_dot2(t, alpha_w, alpha_b, alpha_g): 
  return alpha_w*(alpha_g*beta_w(t, alpha_w, alpha_b, alpha_g) - plague(t)) + 0.001
  
def alphab_dot2(t, alpha_w, alpha_b, alpha_g): 
  return alpha_b*(alpha_g*beta_b(t, alpha_w, alpha_b, alpha_g) - plague(t)) + 0.001

def alphag_dot2(t, alpha_w, alpha_b, alpha_g): 
  return (- Tdot_w2(t, alpha_w, alpha_b, alpha_g) - Tdot_b2(t, alpha_w, alpha_b, alpha_g))
