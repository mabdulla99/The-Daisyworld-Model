############################## defining the Daisyworld model

def Daisyworld(dt, time_array, g0, w0, b0):

alpha_g = np.zeros(len(time_array)) 
alpha_w = np.zeros(len(time_array)) 
alpha_b = np.zeros(len(time_array)) 
temp_p = np.zeros(len(time_array)) 
temp_b = np.zeros(len(time_array))

alpha_g[0] = g0 
alpha_w[0] = w0 
alpha_b[0] = b0

for i in range(0, len(time_array)-1):

  k_1 = dt*alphaw_dot(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i]) 
  k_2 = dt*alphab_dot(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i]) 
  k_3 = dt*alphag_dot(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i])
  
  k_4 = dt*alphaw_dot(time_array[i], alpha_w[i] + (k_1/2), alpha_b[i], alpha_g[i]) 
  k_5 = dt*alphab_dot(time_array[i], alpha_w[i], alpha_b[i] + (k_2/2), alpha_g[i]) 
  k_6 = dt*alphag_dot(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i] + (k_3/2))
  
  k_7 = dt*alphaw_dot(time_array[i], alpha_w[i] + (k_4/2), alpha_b[i], alpha_g[i]) 
  k_8 = dt*alphab_dot(time_array[i], alpha_w[i], alpha_b[i] + (k_5/2), alpha_g[i]) 
  k_9 = dt*alphag_dot(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i] + (k_6/2))
  
  k_10 = dt*alphaw_dot(time_array[i], alpha_w[i] + k_7, alpha_b[i], alpha_g[i]) 
  k_11 = dt*alphab_dot(time_array[i], alpha_w[i], alpha_b[i] + k_8, alpha_g[i]) 
  k_12 = dt*alphag_dot(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i] + k_9)
  
  alpha_w[i+1] = alpha_w[i] + (k_1 + 2*k_4 + 2*k_7 + k_10)/6
  alpha_b[i+1] = alpha_b[i] + (k_2 + 2*k_5 + 2*k_8 + k_11)/6
  alpha_g[i+1] = alpha_g[i] + (k_3 + 2*k_6 + 2*k_9 + k_12)/6
  
  temp_p[i] = ((S*L(time_array[i])* (1-albedo_p(time_array[i], alpha_w[i], alpha_b[i], alpha_g[i])))/(sigma))**(0.25) - 273
  
  temp_b[i] = ((S*L(time_array[i])*(1-0.5))/(sigma))**(0.25) - 273
  
return alpha_w, alpha_b, alpha_g, temp_p, temp_b

a, b, c, d, e = Daisyworld(dt, time, init_g, init_w, init_b)

plt.figure(figsize = (10,8)) 
plt.plot(L(time), a, 'gray', label = 'White Daisies')
7
plt.plot(L(time), b, 'black', label = 'Black Daisies') 
plt.plot(L(time), c, 'brown', label = 'Uncovered Ground') 
plt.title('Covered/Uncovered Land vs. Solar Luminosity', fontsize = 17) 
plt.xlabel('Solar Luminosity', fontsize = 17) 
plt.ylabel('Fraction of World Covered', fontsize = 17) 
plt.legend(fontsize = 15.5)

plt.figure(figsize = (10,8)) 
plt.plot(L(time)[0:-1], d[0:-1], 'red', label = 'Daisyworld') 
plt.plot(L(time)[0:-1], e[0:-1], 'b--', label = 'Barren World') 
plt.title('Mean Planetary Temperature vs. Solar Luminosity', fontsize = 17) 
plt.xlabel('Solar Luminosity', fontsize = 17) 
plt.ylabel('$T_p$ [$^\circ$C]', fontsize = 17) 
plt.legend(fontsize = 15.5) 
