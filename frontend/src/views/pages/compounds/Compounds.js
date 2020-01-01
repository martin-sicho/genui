import React from 'react';
import {useRouteMatch, useParams} from 'react-router-dom';

const Compounds = () => {

  let { path, url } = useRouteMatch();
  let {project} = useParams();

  return (
    <div>
      <span>This is a compounds page. {path}, {url}, {project}</span>
    </div>
  );
};

export default Compounds;
