import numpy as np
import matplotlib.pyplot as plt 


transfer = []
through = []
transfer_peak = []
through_peak = []
s = []
#after file #12, the last items are the db from the sum around, and the ones before are from the max,
#with open("sweep_res13.txt",'r') as f:
with open("sweep_res20.txt",'r') as f:

	#with open("out2.csv",'r') as f:
	for line in f:
		transfer.append(np.float(line.split(',')[-1]))
		through.append(np.float(line.split(',')[-2]))
		transfer_peak.append(np.float(line.split(',')[-3]))
		through_peak.append(np.float(line.split(',')[-4]))
		s.append(np.float(line.split(',')[0]))

transfer = np.array(transfer)
through = np.array(through)
s = np.array(s)

plt.plot(s/1e-6, transfer, label='_nolegend_')
plt.plot(s/1e-6, through, label='_nolegend_')
plt.scatter(s/1e-6, transfer, label='Transfered')
plt.scatter(s/1e-6, through, label='Transmitted')
#plt.vlines([s[np.argmax(transfer)], s[np.argmin(through)]], np.min(transfer), np.max(through))
plt.axvspan((s[np.argmin(through)] - 0.03e-6)/1e-6, (s[np.argmin(through)] + 0.03e-6)/1e-6, facecolor='b', alpha=0.5, label="On State")
plt.axvspan((2e-6 - 0.03e-6)/1e-6, np.max(s)/1e-6, facecolor='g', alpha=0.5, label= "Off State")
plt.legend()
plt.xlabel("$s (\mu m)$")
plt.ylabel("Relative Power (dB)")
plt.ylim(-30, 2)
plt.xlim(np.min(s)/1e-6, np.max(s)/1e-6)
plt.title("Results of Waveguide Separation Sweep")


plt.show()


plt.plot(s, transfer, label='_nolegend_')
plt.plot(s, through, label='_nolegend_')
plt.scatter(s, transfer, label='transfered')
plt.scatter(s, through, label='transmitted')
plt.plot(s, transfer_peak, label='_nolegend_')
plt.plot(s, through_peak, label='_nolegend_')
plt.scatter(s, transfer_peak, label='transfered peak')
plt.scatter(s, through_peak, label='transmitted peak')
plt.legend()
plt.xlabel("$s (\mu m)$")
plt.ylabel("power (dB)")
plt.show()




plt.plot(s, (10*np.ones(len(transfer)))**(transfer/(10*np.ones(len(transfer)))), label='_nolegend_')
plt.plot(s, (10*np.ones(len(through)))**(through/(10*np.ones(len(through)))), label='_nolegend_')
plt.scatter(s, (10*np.ones(len(transfer)))**(transfer/(10*np.ones(len(transfer)))), label='transfered')
plt.scatter(s, (10*np.ones(len(through)))**(through/(10*np.ones(len(through)))), label='transmitted')
#plt.vlines([s[np.argmax(transfer)], s[np.argmin(through)]], np.min(transfer), np.max(through))
plt.legend()
plt.xlabel("$s (\mu m)$")
plt.ylabel("power")
plt.show()


#look at summed power

plt.plot(s, (10*np.ones(len(transfer)))**(transfer/(10*np.ones(len(transfer)))), label='_nolegend_')
plt.plot(s, (10*np.ones(len(through)))**(through/(10*np.ones(len(through)))), label='_nolegend_')
plt.scatter(s, (10*np.ones(len(transfer)))**(transfer/(10*np.ones(len(transfer)))), label='transfered')
plt.scatter(s, (10*np.ones(len(through)))**(through/(10*np.ones(len(through)))), label='transmitted')
plt.scatter(s, (10*np.ones(len(through)))**(through/(10*np.ones(len(through)))) + (10*np.ones(len(transfer)))**(transfer/(10*np.ones(len(transfer)))), label='summed')
#plt.vlines([s[np.argmax(transfer)], s[np.argmin(through)]], np.min(transfer), np.max(through))
plt.legend()
plt.xlabel("$s (\mu m)$")
plt.ylabel("power")
plt.show()







