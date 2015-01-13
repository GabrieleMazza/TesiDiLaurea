#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

using namespace std;

int main()
{
	cout << endl << "CODICI DEI COMUNI DEL VENETO" << endl << endl;
	cout << "Reading sample codes from keys.txt... ";
	
	map<string, int> Keys;

	// Importo i dati
	ifstream file("keys.txt");
	if (file.fail())
	{
		cerr << "ERROR GENERATED!" << endl;
		cerr << "Cannot open keys.txt" << endl;
		system("PAUSE");
		exit(1);
	}

	string str;
	int tmpint;
	string tmpstr;
	int count = 0;
	while (getline(file, str))
	{
		istringstream SSTR(str);
		SSTR >> tmpstr >> tmpint;
		Keys[tmpstr] = tmpint;
		count++;
	}
	file.close();
	cout << "Done (" << count <<" found)." << endl;

	cout << "Detecting codes... ";
	// Importo i dati
	ifstream comuni("in.txt");
	ofstream codes("out.txt");
	if (comuni.fail())
	{
		cerr << "ERROR GENERATED!" << endl;
		cerr << "Cannot open in.txt" << endl;
		system("PAUSE");
		exit(1);
	}
	map<string, int>::iterator it;
	count = 0;
	while (getline(comuni, tmpstr))
	{
		it = Keys.find(tmpstr);
		if (it==Keys.end())
		{
			cerr << "ERROR GENERATED!" << endl;
			cerr << tmpstr << " not defined (row ." <<count+2 << " in excel)" << endl;
			system("PAUSE");
			exit(1);
		}
		else
		{
			codes << Keys[tmpstr] << endl;
		}
		count++;
	}
	
	cout << "Done (" << count << " processed)" << endl;
	system("PAUSE");
	return 0;
}
