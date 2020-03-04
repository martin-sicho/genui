import React from 'react';
import { Button, CardBody, CardFooter } from 'reactstrap';

class ModelFormCardBody extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      formIsSubmitting : false,
    }
  }

  setFormSubmitting = (state) => {
    this.setState({formIsSubmitting : state});
  };

  render() {
    const formIsSubmitting = this.state.formIsSubmitting;
    const Form = this.props.form;

    return (
      <React.Fragment>
        <CardBody className="scrollable">
          <Form
            {...this.props}
            onSubmit={
              (values) => {
                this.setFormSubmitting(true);
                this.props.handleCreate(values);
              }
            }
          />
        </CardBody>
        <CardFooter>
          <Button block form={`${this.props.modelClass}-${this.props.formNameSuffix}-form`} type="submit" color="primary" disabled={formIsSubmitting}>{formIsSubmitting ? "Creating..." : "Create"}</Button>
        </CardFooter>
      </React.Fragment>
    )
  }
}

export default ModelFormCardBody;