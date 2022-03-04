# IC3 format V2

this->data_set->registerVector<double>("X_NO", NO_DATA, false, false);
  this->data_set->registerScalar<int>("NO_FLAG", NO_DATA, false, false);
  this->data_set->registerScalar<int>("FA_FLAG", FA_DATA, false, false);
  this->data_set->registerScalar<int>("CV_FLAG", CV_DATA, false, false);
  this->data_set->registerScalar<int>("CV_PART", CV_DATA, false, false);
  /*global node index, should be conserved during multiple restart process*/
  this->data_set->registerScalar<mint8>("NODE_GB_INDEX", NO_DATA, true, false);
  /*indicates the first two node identified in each cell */
  /* should also be conserved during multiple restart */
  /* This is useful for spectral schemes */
  this->data_set->registerScalar<mint8>("CV_NODE_ROOT_0", CV_DATA, true, false);
  this->data_set->registerScalar<mint8>("CV_NODE_ROOT_1", CV_DATA, true, false);

## restart header

- `UGP_IO_MAGIC_NUMBER` to test byte swap (little endian, big endian)
- version

## section header

order of sections is arbitrary and can be one of the following

- `UGP_IO_NO_FA_CV_NOOFA_COUNTS`: size of mesh{
        msg.str("");
        msg << " Global nno, nfa, ncv: " << header.idata[0] << " "
            << header.idata[1] << " " << header.idata[2];
        Msg(msg.str());
        // ibermejoComment For even larger mesh support, redeclare
        // no_,fa_,cv_count as mint8 and
        // ibermejoComment convert the derived functions using them
        no_count = (int)header.idata[0];
        fa_count = (int)header.idata[1];
        cv_count = (int)header.idata[2];
        noofa_count = header.idata[3];
        msg.str("");
        msg << " noofa_count: " << noofa_count;
        Msg(msg.str());
        this->initRead();
      } break;

- `UGP_IO_NO_CHECK`
- `UGP_IO_FA_CHECK`
- `UGP_IO_CV_CHECK`
- `UGP_IO_NOOFA_I_AND_V`: CSR connectivity of face2node
- `UGP_IO_CVOFA: {`
        readCvofa(header, offset);
- `UGP_IO_FA_ZONE: {`
        readFaZone(header);
- `UGP_IO_CV_PART: {`
        readCvPart(header, offset);
- `UGP_IO_X_NO: {`
        Msg("Reading NO_DATA vector: X_NO ");
        cVector<double>* vec_x_no =
            dynamic_cast<cVector<double>*>(this->data_set->getData("X_NO"));
        this->readNoD3(vec_x_no, header, offset);
        vec_x_no->dumpStats();
- `UGP_IO_DATA:`
- `UGP_IO_I0: {`
        cValue<mint8>* value = dynamic_cast<cValue<mint8>*>(
            this->data_set->getData<cValue<mint8> >(header_name_str));
        if (value != nullptr) {
          msg.str("");
          msg << "reading large int value: " << header_name_str;
          Msg(msg.str());
          value->setValue(header.idata[0]);
          value->setFlag();
          value->dumpRange();
        } else {
          msg.str("");
          msg << "skipping large int value: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_D0: {`
        cValue<double>* value = dynamic_cast<cValue<double>*>(
            this->data_set->getData<cValue<double> >(header_name_str));
        if (value != nullptr) {
          msg.str("");
          msg << " reading double value: " << header_name_str;
          Msg(msg.str());
          value->setValue(header.rdata[0]);
          value->setFlag();
          value->dumpRange();
        } else {
          msg.str("");
          msg << "skipping double value: " << header_name_str << ": "
              << std::endl;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_FA_II1: {`
        cScalar<mint8>* scalar = dynamic_cast<cScalar<mint8>*>(
            this->data_set->getData<cScalar<mint8> >(header_name_str));
        if ((scalar != nullptr) && (scalar->getDatatype() == FA_DATA)) {
          msg.str("");
          msg << "reading mint8 scalar FA_DATA: " << header_name_str;
          Msg(msg.str());
          this->readFaII1(scalar, header, offset);
          scalar->dumpStats();
          scalar->setFlag();

        } else {
          msg.str("");
          msg << "skipping mint8 scalar FA_DATA: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_FA_D1: {`
        cScalar<double>* scalar = dynamic_cast<cScalar<double>*>(
            this->data_set->getData<cScalar<double> >(header_name_str));
        if ((scalar != nullptr) && (scalar->getDatatype() == FA_DATA)) {
          msg.str("");
          msg << "reading double scalar FA_DATA: " << header_name_str;
          Msg(msg.str());
          this->readFaD1(scalar, header, offset);
          scalar->dumpStats();
          scalar->setFlag();

        } else {
          msg.str("");
          msg << "skipping double scalar FA_DATA: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_FA_D3: {`
        cVector<double>* vector = dynamic_cast<cVector<double>*>(
            this->data_set->getData<cVector<double> >(header_name_str));
        if ((vector != nullptr) && (vector->getDatatype() == FA_DATA)) {
          msg.str("");
          msg << "reading double vector FA_DATA: " << header_name_str;
          Msg(msg.str());
          this->readFaD3(vector, header, offset);
          vector->dumpStats();
          vector->setFlag();

        } else {
          msg.str("");
          msg << "skipping double vector FA_DATA: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_CV_II1: {`
        cScalar<mint8>* scalar = dynamic_cast<cScalar<mint8>*>(
            this->data_set->getData<cScalar<mint8> >(header_name_str));
        if ((scalar != nullptr) && (scalar->getDatatype() == CV_DATA)) {
          msg.str("");
          msg << "reading mint8 scalar CV_DATA: " << header_name_str;
          Msg(msg.str());
          this->readCvII1(scalar, header, offset);
          scalar->dumpStats();
          scalar->setFlag();

        } else {
          msg.str("");
          msg << "skipping mint8 scalar CV_DATA: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_CV_D1: {`
        cScalar<double>* scalar = dynamic_cast<cScalar<double>*>(
            this->data_set->getData<cScalar<double> >(header_name_str));
        if ((scalar != nullptr) && (scalar->getDatatype() == CV_DATA)) {
          msg.str("");
            msg << "reading double scalar CV_DATA: " << header_name_str;
          Msg(msg.str());
          this->readCvD1(scalar, header, offset);
          scalar->dumpStats();
          scalar->setFlag();
        } else {
          msg.str("");
          msg << "skipping double scalar CV_DATA: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_CV_D3: {`
        cVector<double>* vector = dynamic_cast<cVector<double>*>(
            this->data_set->getData<cVector<double> >(header_name_str));
        if ((vector != nullptr) && (vector->getDatatype() == CV_DATA)) {
          msg.str("");
          msg << "reading double vector CV_DATA: " << header_name_str;
          Msg(msg.str());
          this->readCvD3(vector, header, offset);
          vector->dumpStats();
          vector->setFlag();

        } else {
          msg.str("");
          msg << "skipping double vector CV_DATA: " << header_name_str;
          Msg(msg.str());
        }
      } break;
- `UGP_IO_CV_D33`: cell based and double type 3x3 tensor data
- `UGP_IO_NO_II1`: node based and int scalar data
- `UGP_IO_NO_D1`: node based and double scalar data
- `UGP_IO_NO_D3`: node based and double vector data
- `UGP_IO_EOF`: end of file

## 

## mesh data

order of sections is arbitrary

### `NOOFA` face2node connectivity

- header with number of nodes, faces and cells