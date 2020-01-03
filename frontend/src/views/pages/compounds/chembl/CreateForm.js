import { Field, FieldArray, Formik } from 'formik';
import * as Yup from 'yup';
import {
  Button,
  CardBody,
  CardFooter, Col,
  Form,
  FormGroup,
  Input,
  Label, Row,
  UncontrolledAlert,
} from 'reactstrap';
import React from 'react';

const ErrorMessage = ({ name }) => (
  <Field
    name={name}
  >
    {
      ({field, form, meta}) => {
        const error = meta.error;
        const touch = meta.touched;
        if (touch && error) {
          return <UncontrolledAlert color="danger">{error}</UncontrolledAlert>;
        }

        return null;
      }
    }
  </Field>
);

export function ChEMBLCreateForm(props) {

  return (
    <React.Fragment>
      <CardBody className="scrollable">
        <Formik
          initialValues={{
            name: 'New Compound Set Name',
            description: 'Detailed description of the compound set...',
            targets : []
          }}
          validationSchema={Yup.object(
            {
              name: Yup.string()
                .max(256, 'Name must be less than 256 character long.')
                .required('Name is required.'),
              description: Yup.string()
              .max(10000, 'Description must be 10,000 characters or less'),
              targets: Yup.array().of(Yup.string().required("Empty target not allowed")).required("Must provide at least one target.")
            })
          }
          onSubmit={ values => props.handleCreate(values)}
        >
          {
            formik => (
              <Form id="chembl-download-form" onSubmit={formik.handleSubmit} className="unDraggable">
                <FormGroup>
                  <Label htmlFor="name">Name</Label>
                  <Field name="name" as={Input} type="text"/>
                </FormGroup>
                <ErrorMessage name="name"/>
                <FormGroup>
                  <Label htmlFor="description">Description</Label>
                  <Field name="description" as={Input} type="textarea"/>
                </FormGroup>
                <ErrorMessage name="description"/>
                  <FieldArray
                    name="targets"
                    render={(arrayHelpers) => (
                      <React.Fragment>
                      <FormGroup>
                        <Label htmlFor="Targets">Targets</Label>
                        {formik.values.targets && formik.values.targets.length > 0 ? (
                          formik.values.targets.map((target, index) => (
                            <React.Fragment key={index}>
                              <FormGroup>
                                <Row>
                                  <Col sm={10}>
                                    <Field name={`targets.${index}`} as={Input} type="text" />
                                  </Col>
                                  <Col sm={2}>
                                    <Button
                                      type="button"
                                      onClick={() => arrayHelpers.remove(index)} // remove a friend from the list
                                    >
                                      -
                                    </Button>
                                    <Button
                                      type="button"
                                      onClick={() => {arrayHelpers.insert(index, "")}} // insert an empty string at a position
                                    >
                                      +
                                    </Button>
                                  </Col>
                                </Row>
                              </FormGroup>
                              <ErrorMessage name={`targets.${index}`}/>
                            </React.Fragment>
                          ))
                        ) : (
                          <Button type="button" onClick={() => arrayHelpers.push("")}>
                            {/* show this when user has removed all friends from the list */}
                            Add a ChEMBL Target
                          </Button>
                        )}
                      </FormGroup>
                    {typeof formik.errors.targets === 'string' ?
                      <UncontrolledAlert color="danger">{formik.errors.targets}</UncontrolledAlert>
                      : null}
                      </React.Fragment>
                    )}
                  />
              </Form>
            )
          }
        </Formik>
      </CardBody>
      <CardFooter>
        <Button form="chembl-download-form" type="submit" color="primary">Create</Button>
      </CardFooter>
    </React.Fragment>
  );
}