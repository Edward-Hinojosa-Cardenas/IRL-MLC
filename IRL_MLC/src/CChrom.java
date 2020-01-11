public class CChrom {
	public double[] rule;
	public CChrom(int sizeRule) {
		// TODO Auto-generated constructor stub
		rule = new double[sizeRule];
	}
	public CChrom(double new_rule[]) {
		// TODO Auto-generated constructor stub
		rule = new double[new_rule.length];
		for(int i = 0; i < new_rule.length; i++) {
			rule[i] = new_rule[i];
		}
	}
	public int getSize() {
		return rule.length;
	}
	public double getValu(int index) {
		return rule[index];
	}
}
