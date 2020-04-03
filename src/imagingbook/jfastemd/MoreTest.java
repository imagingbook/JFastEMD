package imagingbook.jfastemd;

public class MoreTest {

	public static void main(String[] args) {
		for (int i = 0; i < getLimit(); i++) {
			System.out.println(i);
		}

	}
	
	static int getLimit() {
		System.out.println("get");
		return 10;
	}

}
