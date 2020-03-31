import React from "react";
import withUnmounted from '@ishawnwang/withunmounted';

class ComponentWithObjects extends React.Component {
  abort = new AbortController();
  hasUnmounted = false;

  constructor(props) {
    super(props);

    this.objectListRoot = this.props.objectListURL;
    this.emptyClassProperty = props.emptyClassName;

    if (!this.emptyClassProperty) {
      throw new Error("Unspecified empty class name for ComponentWithObjects. This is required. Specify a default name if class cannot be determined from data.");
    }

    const objects = {};
    objects[this.emptyClassProperty] = [];
    this.state = {
      fetchUpdates : false,
      objects : objects,
    }
  }

  componentDidMount() {
    this.setState({fetchUpdates : true})
  }

  componentWillUnmount() {
    this.abort.abort();
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.currentProject && this.state.fetchUpdates) {
      this.fetchUpdates(this.props.currentProject);
    }
  }

  fetchUpdates = (project) => {
    const params = new URLSearchParams();
    params.append('project_id', project.id);
    fetch(this.objectListRoot.toString() + "?" + params.toString(), {signal : this.abort.signal, credentials: "include",})
      .then(response => response.json())
      .then(this.getObjects)
      .catch(
        (error) => console.log(error)
      );
  };

  getObjects = (data) => {
    const objects = {};
    objects[this.emptyClassProperty] = [];
    for (const obj of data) {
      if (!obj.className) {
        obj.className = this.emptyClassProperty
      }
      if (!objects.hasOwnProperty(obj.className)) {
        objects[obj.className] = [];
      }
      objects[obj.className].push(obj);
    }

    if (this.hasUnmounted) {
      return
    }
    this.setState({
      objects : objects,
      fetchUpdates : false
    })
  };

  handleAddObject = (className, data) => {
    this.setState(prevState => {
      if (!prevState.objects.hasOwnProperty(className)) {
        prevState.objects[className] = []
      }
      prevState.objects[className].push(data);
      return {
        objects : prevState.objects
      };
    });
  };

  handleAddObjectList = (className, objectList, overwrite=false) => {
    this.setState((prevState) => {
      const old_objects = prevState.objects;
      if (old_objects.hasOwnProperty(className)) {
        if (overwrite) {
          old_objects[className] = objectList;
        } else {
          old_objects[className] = old_objects[className].concat(objectList);
        }
      } else {
        old_objects[className] = objectList;
      }
      return {
        objects : old_objects
      }
    });
  };

  deleteFromState = (className, object) => {
    this.setState(prevState => {
      const objects = prevState.objects[className];
      const idx_del = objects.findIndex(item => item.id === object.id);
      objects.splice(idx_del, 1);
      if (objects.length === 0) {
        prevState.objects[className] = []
      }
      return {
        objects : prevState.objects
      };
    });
  };

  handleObjectDelete = (className, object) => {
    if (this.props.customDelete) {
      this.props.customDelete(className, object);
      this.deleteFromState(className, object);
    } else {
      fetch(this.objectListRoot.toString() + object.id + '/', {method: 'DELETE', credentials: "include",})
        .then(
          () => {
            this.deleteFromState(className, object)
          }
        ).catch(
        (error) => console.log(error)
      )
      ;
    }
  };

  requestObjectUpdate = () => {
    this.setState({fetchUpdates: true})
  };

  render() {
    if (this.state.objects === null) {
      return <div>Loading...</div>
    }

    return (
      <React.Fragment>
        {this.props.render(this.state.objects, this.handleAddObjectList, this.handleAddObject, this.handleObjectDelete, this.requestObjectUpdate)}
      </React.Fragment>
    )
  }
}

export default withUnmounted(ComponentWithObjects);