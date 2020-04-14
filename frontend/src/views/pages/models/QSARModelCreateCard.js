import * as Yup from 'yup';
import { MolsetActivitiesSummary, ModelCardNew, SimpleDropDownToggle } from '../../../genui';
import React from 'react';
import { QSARExtraFields, QSARTrainingFields, QSARValidationFields } from './QSARModelFormFields';
import { Button, CardBody, CardHeader, Col, Row, CardFooter } from 'reactstrap';

function EndpointSelector(props) {
  return <MolsetActivitiesSummary {...props} selectable={true} message="Choose the desired activity endpoint by clicking the corresponding row in the table below. The chosen activity type from the given activity set will be used as the output variable for the resulting model."/>
}

export default function QSARModelCreateCard (props) {
  let molsets = [];
  Object.keys(props.compoundSets).forEach(
    (key) => molsets = molsets.concat(props.compoundSets[key])
  );

  const [molset, setMolset] = React.useState(null);
  const [endpointData, setEndpointData] = React.useState(null);
  const [dataReady, setDataReady] = React.useState(false);

  const trainingStrategyInit = {
    activityThreshold : 6.5,
    descriptors: [props.descriptors[0].id],
  };
  const validationStrategyInit = {
    cvFolds: 10,
    validSetSize: 0.2,
  };
  const extraParamInit = {
    molset: molset ? molset.id : undefined,
    predictionsType: endpointData ? endpointData.type.value : undefined,
    predictionsUnits: ""
  };

  const trainingStrategySchema = {
    activityThreshold: Yup.number().min(0, 'Activity threshold must be zero or positive.').required('Activity threshold is a required parameter.'),
    descriptors: Yup.array().of(Yup.number().positive('Descriptor set ID must be a positive integer.')).required('You need to supply one or more descriptor sets for training.'),
  };
  const validationStrategySchema = {
    cvFolds: Yup.number().integer().min(0, 'Number of CV folds must be at least 0.'),
    validSetSize: Yup.number().min(0.0, 'Validation set size must be at least 0.0.').max(1.0,'Validation set size is expressed as a fraction, which needs to be less than 1.0.'),
  };

  const extraParamsSchema = {
    molset: Yup.number().integer().positive('Molecule set ID must be a positive integer.').required('You need to supply a training set of compounds.'),
    predictionsType: Yup.string().required('Predictions activity type has to be set.').min(1, "Predictions type cannot be empty.").max(128, 'Predictions type name cannot be longer than 128 characters.'),
    predictionsUnits: Yup.string().max(128, 'Predictions units name cannot be longer than 128 characters.')
  };

  if (endpointData) {
    trainingStrategyInit.activityType = endpointData.type.id;
    trainingStrategyInit.activitySet = endpointData.activitySet.id;

    trainingStrategySchema.activityType = Yup.number().integer().positive('Activity type ID must be a positive integer.').required('You need to select an activity type for modelling.');
    trainingStrategySchema.activitySet = Yup.number().integer().positive('Activity set ID must be a positive integer.').required('You need to supply a set of activities to use for modelling.');
  }

  return !dataReady ? (
    <React.Fragment>
      <CardHeader>QSAR Training Set and Activity Endpoint</CardHeader>
      <CardBody className="scrollable">
        <SimpleDropDownToggle
          items={molsets}
          onSelect={setMolset}
          message={() => <p>You have to choose a training set. You can choose any compound set in the current project.</p>}
          title="Choose Training Set"
          header="Available Compound Sets"
        />

        {
          molset ? (
            <React.Fragment>
              <hr/>
              <Row>
                <Col sm={12}>
                  <EndpointSelector {...props} molset={molset} onSelect={setEndpointData}/>
                </Col>
              </Row>
            </React.Fragment>
          ) : null
        }
      </CardBody>

      {
        molset && endpointData ? (
          <CardFooter>
            <Row>
              <Col sm={8}>
                <p>Selected endpoint: {endpointData.name}</p>
              </Col>
              <Col sm={4}>
                <Button block color="primary" onClick={() => setDataReady(true)}>Choose Model Parameters</Button>
              </Col>
            </Row>
          </CardFooter>
        ) : null
      }
    </React.Fragment>
  ) : (
    <ModelCardNew
      {...props}
      molsets={[molset]}
      activitySets={[endpointData.activitySet]}
      activityTypes={[endpointData.type]}
      endpointData={endpointData}
      trainingStrategyInit={trainingStrategyInit}
      validationStrategyInit={validationStrategyInit}
      extraParamsInit={extraParamInit}
      trainingStrategySchema={trainingStrategySchema}
      validationStrategySchema={validationStrategySchema}
      extraParamsSchema={extraParamsSchema}
      trainingStrategyFields={QSARTrainingFields}
      validationStrategyFields={QSARValidationFields}
      extraFields={QSARExtraFields}
      onValuesInit={(values, state) => {
        if (state.modes) {
          const mode = state.modes[0];
          if (mode.name === "classification") {
            values.predictionsType = "Active Probability"
          }
        }
        return values;
      }}
      prePost={(data) => {
        if (data.predictionsUnits === "") {
          data.predictionsUnits = null;
        }
        return data;
      }}
    />
  )
}