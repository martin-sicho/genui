import React from 'react';
import * as Yup from 'yup';
import { SimpleDropDownToggle } from '../../index';
import { CardBody } from 'reactstrap';

class ModelFormRenderer extends React.Component {
  CTYPE_TO_VALIDATOR = {
    string: Yup.string().required(),
    integer: Yup.number().required(),
    float: Yup.number().required(),
    bool: Yup.bool().required()
  };

  constructor(props) {
    super(props);

    this.chosenAlgorithm = this.props.chosenAlgorithm;
    this.parameters = !this.props.omitAlgParams ? this.props.chosenAlgorithm.parameters : null;
    this.enableFileUploads = this.props.enableFileUploads;
    this.disabledModelFormFields = this.props.disabledModelFormFields ? this.props.disabledModelFormFields : [];

    this.state = {
      metrics: !this.props.omitValidation && this.props.metrics ? this.initMetrics(this.props.metrics) : null,
      modes: this.props.chosenAlgorithm.validModes,
      initialValues: null,
      schema: null,
      formDataReady: false,
    }
  }

  componentDidMount() {
    if (this.state.modes.length === 1) {
      this.initFormData();
    }
  }

  initMetrics = (metrics) => {
    const ret = [];
    metrics.forEach(metric => {
      if (metric.validAlgorithms.length === 0 || metric.validAlgorithms.find(alg => this.chosenAlgorithm.id === alg)) {
        ret.push(metric);
      }
    });
    return ret;
  };

  initFormData = () => {
    this.setState({
      initialValues: this.generateInit(),
      schema: this.generateSchema(),
      formDataReady: true,
    })
  };

  handleModeSelect = (mode) => {
    if (this.state.metrics) {
      const metrics = [];
      this.state.metrics.forEach(metric => {
        if (metric.validModes.find(item => mode.id === item.id)) {
          metrics.push(metric);
        }
      });
      this.setState({
        metrics: metrics,
        modes: [mode]
      }, this.initFormData)
    } else {
      this.setState({
          modes: [mode]
        }, this.initFormData
      )
    }
  };

  generateInit = () => {
    const trainingStrategyDefaultInit = {
      algorithm: this.chosenAlgorithm.id,
      mode: this.state.modes.length > 0 ? this.state.modes[0].id : [],
    };
    const trainingStrategyInit = Object.assign(trainingStrategyDefaultInit, this.props.trainingStrategyInit ? this.props.trainingStrategyInit : {});
    const validationStrategyDefaultInit = this.state.metrics && !this.disabledModelFormFields.includes('validationStrategy.metrics') ? {
      metrics: this.state.metrics.length > 0 ? [this.state.metrics[0].id] : []
    } : {};
    const validationStrategyInit = Object.assign(validationStrategyDefaultInit, this.props.validationStrategyInit ? this.props.validationStrategyInit : {});

    let initialValues = {
      name: `New ${this.chosenAlgorithm.name} Model`,
      description: '',
      project: this.props.project.id,
      trainingStrategy: trainingStrategyInit
    };
    if (this.enableFileUploads) {
      initialValues.modelFile = undefined;
    }

    // validation
    if (Object.keys(validationStrategyInit).length !== 0 && validationStrategyInit.constructor === Object) {
      initialValues.validationStrategy = validationStrategyInit;
    }

    // extra parameters
    initialValues = Object.assign(initialValues, this.props.extraParamsInit ? this.props.extraParamsInit : {});

    // default parameters
    if (this.parameters) {
      const parameterDefaults = {};
      for (const param of this.parameters) {
        parameterDefaults[param.name] = param.defaultValue;
      }
      initialValues.trainingStrategy.parameters = parameterDefaults;
    }

    if (this.props.onValuesInit) {
      initialValues = this.props.onValuesInit(initialValues, this.state);
    }

    // console.log(initialValues);

    return initialValues;
  };

  generateSchema = () => {
    const trainingStrategyDefault = {
      algorithm: Yup.number().integer().positive("Algorithm ID needs to be a positive number").required('Algorithm ID must be supplied'),
      mode: Yup.number().integer()
        .max(256, 'Mode must be 256 characters or less.').required('You must specify a mode.'),
    };

    // parameters
    if (this.parameters) {
      const parameterValidators = {};
      for (const param of this.parameters) {
        parameterValidators[param.name] = this.CTYPE_TO_VALIDATOR[param.contentType]
      }
      trainingStrategyDefault.parameters = Yup.object().shape(parameterValidators);
    }

    const trainingStrategy = Object.assign(trainingStrategyDefault, this.props.trainingStrategySchema);
    // the main schema object
    let validationObj = {
      name: Yup.string()
        .max(256, 'Name must be less than 256 characters long.')
        .required('Name is required.'),
      description: Yup.string()
        .max(10000, 'Description must be 10,000 characters or less.'),
      project: Yup.number().integer().positive("Project ID needs to be a positive number").required('Project ID must be supplied'),
      trainingStrategy: Yup.object().shape(
        trainingStrategy
      )
    };
    if (this.enableFileUploads) {
      validationObj.modelFile = Yup.mixed().required("Model file is required.");
    }

    // validation added only if there is something to add
    const validationStrategyDefault = this.state.metrics && !this.disabledModelFormFields.includes('validationStrategy.metrics') ? {
      metrics: Yup.array().of(Yup.number().positive('Metric ID must be a positive integer.')).required('You need to supply at least one metric for validation.'),
    } : {};
    const validationStrategy = Object.assign(validationStrategyDefault, this.props.validationStrategySchema);
    if (Object.keys(validationStrategy).length !== 0 && validationStrategy.constructor === Object) {
      validationObj.validationStrategy = Yup.object().shape(
        validationStrategy
      );
    }

    // extra parameters
    validationObj = Object.assign(validationObj, this.props.extraParamsSchema);

    if (this.props.onSchemaInit) {
      validationObj = this.props.onSchemaInit(validationObj, this.state);
    }

    // console.log(validationObj);
    return Yup.object().shape(validationObj);
  };

  render() {
    if (!this.state.formDataReady) {
      return (
        <CardBody>
          <SimpleDropDownToggle
            items={this.state.modes}
            onSelect={this.handleModeSelect}
            message={() => (
              <p>
                {this.chosenAlgorithm.name} model can be built in {this.state.modes.length} modes. {!this.enableFileUploads ? "Select the desired one below." : "Select the correct mode for the uploaded model below."}
              </p>
            )}
            title="Choose Mode"
            header="Available Modes"
          />
        </CardBody>
      )
    } else {
      const FormComponent = this.props.component;
      return (
        <FormComponent
          {...this.props}
          initialValues={this.state.initialValues}
          validationSchema={this.state.schema}
          modes={this.state.modes}
          metrics={this.state.metrics}
          parameters={this.parameters}
        />
      );
    }
  }
}

export default ModelFormRenderer;