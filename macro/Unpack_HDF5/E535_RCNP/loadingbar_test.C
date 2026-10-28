#include <stdio.h>
#include <windows.h>
using namespace std;

void loadingBar()
{
	system("color0A");

	char a = 177, b = 219;

	printf("\n\n\n\n");
	printf("\n\n\n\n\t\t\t\t\t" + "Loading...\n\n");
	printf("\t\t\t\t\t");

	for (int i = 0; i< 26; i++)
		printf("%c", a);

	printf("\r");
	printf("\t\t\t\t\t");

	for (int i = 0; i < 26; i++) {
		printf("%c", b);

		Sleep(1000);
	}
}

int main()
{
	loadingBar();
	return 0;
}
