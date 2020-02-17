import React from 'react';
import { CardHeader } from 'reactstrap';
import ModelFormRenderer from './ModelFormRenderer';
import ModelFormCardBody from './ModelFormCardBody';
import FormikModelForm from './FormikModelForm';

class ModelCardNew extends React.Component {

  newModelFromFormData = (data) => {
    fetch(
      this.props.listURL
      , {
        method: 'POST'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        }
      }
    ).then((data) => this.props.handleResponseErrors(data, "Request failed. Data wrong or incomplete?")).then(
      data => {
        this.props.handleAddModel(this.props.modelClass, data)
      }
    ).catch(
      error => console.log(error)
    );
  };

  render() {
    return (
      <React.Fragment>
        <CardHeader>Create New {this.props.chosenAlgorithm.name} Model</CardHeader>
        <ModelFormRenderer
          {...this.props} // all these props will be passed down to the component
          component={props => <ModelFormCardBody {...props} form={FormikModelForm}/>} // this is what should draw the formik form and pass the renderer props to it
          handleCreate={this.newModelFromFormData} // this is the method used to process the parsed data from the form
          project={this.props.currentProject} // this is required
          formNameSuffix="create" // this is required
        />
      </React.Fragment>
    )
  }
}

export default ModelCardNew;