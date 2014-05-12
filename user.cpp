//this file is the interface between your feature class and crfparser
//the following varibles can be used directly
//_fpool->_cols: column number of traning file
//base_feature::_labels: dependency relation type list
#include <iostream>
#include "crfparser.h"
/////////////////////////////////////////
// user added code
/////////////////////////////////////////
#include "templet_feature.h"
/////////////////////////////////////////
// end
/////////////////////////////////////////
using namespace std;
bool crfparser::learn_api(char *training_file){
//feature 1: templet feature
/////////////////////////////////////////
// user added code
/////////////////////////////////////////
	templet_feature *tf=new templet_feature(1,base_feature::_labels.size());
	tf->load_templet("template");
	tf->generate_feature_candidate(training_file);
	tf->write_model("tf");
	delete tf;
	tf=new templet_feature(1,base_feature::_labels.size());
	tf->load_model("tf",false,true);
	_fpool->push_back(tf);
	cout<<"feature number: "<<_fpool->_feature_num<<endl;
/////////////////////////////////////////
// end
/////////////////////////////////////////
	_fpool->generate_feature(training_file);
	_fpool->pop();


/*
//feature x:
/////////////////////////////////////////
// user added code
/////////////////////////////////////////
	_feature_x *xf ...
	...
	_fpool->push_back(_feature_x);
/////////////////////////////////////////
// end
/////////////////////////////////////////
	cout<<"feature number: "<<_fpool->_feature_num<<endl;
	_fpool->generate_feature(training_file);
	_fpool->pop();
*/

	return true;
}


bool crfparser::model_api(bool **fmap,char *training_file, int action){
	templet_feature *tf=new templet_feature(1,base_feature::_labels.size());
	if(action==0){//refine, save, then delete
/////////////////////////////////////////
// user added code
/////////////////////////////////////////
//feature 1: templet feature
		if(!tf->load_model("tf",true,false))
			return false;
		tf->refine_feature(fmap[0]);
		tf->write_model("tf");
		_fpool->pop();
/////////////////////////////////////////
// end
/////////////////////////////////////////
	}else if(action==1){//generate feature, and delete
/////////////////////////////////////////
// user added code
/////////////////////////////////////////
//feature 1: templet feature
		if(!tf->load_model("tf",false,true))
			return false;
		_fpool->push_back(tf);
/////////////////////////////////////////
// end
/////////////////////////////////////////
		_fpool->generate_feature(training_file);
		_fpool->pop();
	}else if(action==2){//just load
/////////////////////////////////////////
// user added code
/////////////////////////////////////////
//feature 1: templet feature
		if(!tf->load_model("tf",false,true))
			return false;
		_fpool->push_back(tf);
/////////////////////////////////////////
// end
/////////////////////////////////////////
	}
	return true;
}
